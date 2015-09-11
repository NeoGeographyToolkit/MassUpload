#!/usr/bin/env python
import httplib
import hashlib
import os
import sys
import time
import shlex
import subprocess
import argparse
import multiprocessing, logging
import Queue
import tempfile
import itertools

''' This is a streamlined version of the GoogleNasa upload tool.  It will recursively upload
    a folder to a Google Cloud Bucket location.  Repeat uploads are limited to files modified
    after a certain timestamp.
'''

# TODO: The multi-threading seems to have permission problems

logger = logging.getLogger()
logger.setLevel(logging.INFO)


def gsutil(*args):
    args = (gsutil_path,) + args
    print " ".join(args)
    return subprocess.check_call(args)

def md5_digest(filepath):
    md5 = hashlib.md5()
    with open(filepath, 'rb') as f:
        for bytes in iter(lambda: f.read(128*md5.block_size), ''):
            md5.update(bytes)
    return md5.hexdigest()

def touch(fname, times=None, time=None):
    if time and not times:
        # set atime and mtime both equal to time
        times = (time, time)
    with file(fname, 'a'):
        os.utime(fname, times) # times is None for current time, or a tuple of timestamps: (atime, mtime)

def list_modified_files(root_path, reference_file):
    """
    Search the filesystem, starting at root_path.
    Return a list of all files in the subtree with modification date
    greater than that of a reference file
    """
    if os.path.exists(reference_file):
        logger.info( "Searching for files modified since: %s" % time.ctime(os.path.getmtime( reference_file )) )
        args = ('find', root_path, '-newer', reference_file)
    else:
        logger.info( "No mtime reference file available, finding all files in folder.")
        args = ('find', root_path)

    logging.debug(' '.join(args) )
    args = shlex.split( ' '.join( args ) )
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    for line in p.stdout.readlines():
        path = line.strip()
        if not os.path.isdir(path):
            relpath = os.path.relpath(path, root_path)
            #logger.info("Found modified file: %s" % relpath)
            yield relpath

def show_modified(options):
    sync_timestamp_file = options.sync_timestamp_file
    if not os.path.exists( sync_timestamp_file ):
        print "%s does not exist.  Aborting." % sync_timestamp_file
        sys.exit(1)

    for file in list_modified_files(options.source_dir, sync_timestamp_file):
        print file
        sys.stdout.flush()


def sync_file(path, bucket, start_path, options):
    logger.info( "SYNC " + path + "...")
    if start_path not in path: # start_path may be a base path or an absolute path, make it absolute.
        path = os.path.join(start_path, path)
    relpath = os.path.relpath(path, start_path)
    
    if options.prepend_path:
        relpath = os.path.join(options.prepend_path, relpath)

    logger.debug("UPLOADING " + path)
    gsutil('-m', 'cp', path, "gs://"+bucket+"/"+relpath)
    #gsutil('cp', path, "gs://"+bucket+"/"+relpath)
    logger.debug("UPLOADED " +path)

def transfer_chunk(paths_in, bucket, root_path, options):
    """  
    Copy a given sequence of filenames to a GCS bucket such that the object keys are all prefixed by 
    each file's path relative to root_path.  gsutil can't do this on it's own so we need to 
    create a temporary subtree on the filesystem representing the relative path tree and 
    then recursively transfer that subtree.
    """
    tmpdir = tempfile.mkdtemp()

    # Filter out files that don't exist locally
    paths = []
    for p in paths_in:
        fullPath = os.path.join(root_path, p)
        if os.path.exists(fullPath):
            paths.append(p)

    # Assuming they're either all absolute or all relative
    relpaths = []
    if root_path in paths[0]: # Compute relative paths from absolute
        for path in paths:
            relpath = os.path.relpath(path, root_path)
            relpaths.append(relpath)
    else: # Already relative
        relpaths = paths

    logger.info("Transferring: %s and %d more..." % (relpaths[0], len(paths)-1) )


    # Make all the sub-folders in the temp directory
    for relpath in relpaths:
        try:
            if options.prepend_path: # Handle prepend option
                relpath = os.path.join(options.prepend_path, relpath)
            temp_folder = os.path.join(tmpdir, os.path.dirname(relpath) )
            os.makedirs(temp_folder)
        except OSError: pass

    # Fill up the temporary directory with symlinks to all the real files
    #links = []
    for relpath in relpaths:
        if os.path.isdir(os.path.join(root_path, relpath)): # Skip folders, they don't need links.
            continue
        original_path = os.path.join(root_path, relpath)
        if options.prepend_path: # Handle prepend option
            relpath = os.path.join(options.prepend_path, relpath)
        symlink_target = os.path.join(tmpdir, relpath)
        # some failure modes have resulted in the link already existing... duplicate records? reissued tmp dirs?
        if os.path.exists(symlink_target):
            os.unlink(symlink_target) 
        os.symlink( original_path, symlink_target )
        #links.append( symlink_target )

    gsutil('-m', 'cp', '-R', os.path.join(tmpdir, '*'), 'gs://'+bucket)

    # Clean up
    subprocess.check_call( ('rm', '-rf', tmpdir) )


def worker(task_q, finished, options):
    '''First worker thread function, this one uploads files to google.
       The output of this thread is sent straight to db_reporter threads.'''
    while True:
        if finished.is_set():
            break
        try:
            # Note: Path may contain many individual paths!
            path, dest_bucket, start_path = task_q.get(False, 0.1) # retry after a timeout
            logger.debug("unqueued (%s, %s, %s)" % (path, dest_bucket, start_path) )
        except Queue.Empty:
            continue
        try:
            if options.by_chunk:
                transfer_chunk(path, dest_bucket, start_path, options)
            else:
                sync_file(path, dest_bucket, start_path, options)
        except Exception, e:
            # Log exceptions that occur and flag this task as failed.
            print str(e)
            logger.info("Caught exception uploading: " + str(path) + '\n')
            raise Exception('FAIL!!!')
        task_q.task_done()


def sync_parallel(source_path_iterator, options):
    '''Note that source_path_iterator can yield files or lists of files!'''
    global connection
    # make uploads world-readable
    gsutil( "defacl", "set",  "public-read", "gs://"+options.dest_bucket) 
    response = None

    # Set up job queues
    task_q   = multiprocessing.JoinableQueue(options.max_queue_size)
    finished = multiprocessing.Event()
    # Create file upload worker threads
    print 'Starting ' + str(options.num_processes) + ' worker processes.'
    subprocesses = [ multiprocessing.Process(target=worker, args=(task_q, finished, options)) for i in range(options.num_processes) ]
    # Start up all the jobs
    for p in subprocesses:
        p.start()
    for paths in source_path_iterator:
        task_q.put( (paths, options.dest_bucket, options.source_dir) )
    # Wait for the jobs to complete.
    task_q.join()
    finished.set()
    logger.info("Finished.  Waiting for worker threads to join")
    for p in subprocesses:
        p.join()
    logger.info("Finished cleaning up worker threads.")

def get_chunks(path_generator):
    """
    Group the output of path_generator into tuples of length options.chunk_size.
    """
    while True:
        paths = []
        for i in range( options.chunk_size):
            try:
                nextVal = path_generator.next()
                # This should only get individual strings as input!
                if not isinstance(nextVal, basestring):
                    print nextVal
                    print len(nextVal)
                    raise Exception('Not a string!')
                paths.append( nextVal )
            except StopIteration: 
                break
        if len(paths) < 1: 
            break
        yield tuple(paths)

def search_for_gsutil():
    # Search for gsutil executable in the $PATH
    for path in os.environ['PATH'].split(':'):
        testpath = os.path.join(path, 'gsutil')
        if os.path.exists(testpath):
            if os.access(testpath, os.X_OK):  # test executability
                return testpath
                break
            else:
                descend = os.path.join(testpath, 'gsutil')
                if os.path.exists(descend) and os.access(descend, os.X_OK):
                    return descend
                    break
    else:
        raise Exception("Couldn't find gsutil in your $PATH.  Please add it or specify it with the --gsutil-path option.")
    
def ensure_exists(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['sync',  'sync-parallel', 'show-modified'])

    parser.add_argument('--dir',         dest='source_dir', metavar='root source directory', 
                                           default='/byss/docroot/smcmich1/hrscMosaicKml/')
    parser.add_argument('--bucket',      dest='dest_bucket', default='hrsc_map_storage')
    parser.add_argument('--prepend-path', dest='prepend_path', default='',
                                          help='Prepend this path to the uploaded files')
    parser.add_argument('--gsutil-path', dest='gsutil_path')
    parser.add_argument('--chunk-size',    dest='chunk_size', type=int, default=1,
         help="For parallel sync, optimize by passing entire chunks to gsutil.")
    parser.add_argument('-p', type=int,  dest='num_processes', 
         help='Number of subprocesses to spawn for sync-parallel', default=8)
    parser.add_argument('--max_queue_size', type=int, help="Size of the task queue", default=20)
    parser.add_argument('--debug',       action='store_true', default=False, help="Turn on debug logging.")

    options = parser.parse_args()
    print "%s --> %s" % (options.source_dir, options.dest_bucket)

    # derived options
    options.by_chunk = options.chunk_size > 1
    options.inventory_path      = options.source_dir
    options.sync_timestamp_file = os.path.join(options.inventory_path, 'last_sync' ) # this file exists solely to mark the time the last sync completed.
    
    if options.debug:
        logger.setLevel(logging.DEBUG)
        multiprocessing.log_to_stderr().setLevel(multiprocessing.SUBDEBUG)


    # look for gsutil in some reasonable places if it's not provided at the command line
    global gsutil_path
    if options.gsutil_path:
        gsutil_path = options.gsutil_path
    else:
        gsutil_path = search_for_gsutil()

    start_time = time.time()

    logging.debug("Start mode switching.")

    if options.command == 'show-modified':
        show_modified(options)
    elif options.command == 'sync-parallel' or options.command == 'sync':
        path_sources = []

        logging.debug("use_modtime: Checking for modified files.")
        # Make an iterator to all the files which have been modified
        modified_files_iter = list_modified_files(options.source_dir, options.sync_timestamp_file)
        path_sources.append(modified_files_iter)

        # Merge all the path sources iterators into a single iterator
        path_generator = itertools.chain( *path_sources )

        if options.by_chunk:
            path_generator = get_chunks( path_generator )
            
        sync_parallel(path_generator, options)

        #if options.use_modtime: # only touch the timestamp file if we used modtime for the sync
        touch( options.sync_timestamp_file ) # update last sync timestamp

    else:
        assert False # argparse shouldn't let us get to this condition

    print '---=== deploy.py script is finished!  ===---'
