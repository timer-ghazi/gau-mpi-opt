#!/usr/bin/env -S python -m mpi4py
# 
#
# --- Embarrassingly Parallel Gaussian Geometry Optimization with an External Code --- 
#                    High-Level Energies + Gaussian Optimizer 
#                        Q.T., Marquette U., 2018--2024
#
#
import sys
import re   
import os
import os.path
import molecule
from molecule import Molecule,Atom
import time
import signal
import subprocess
from mpi4py import MPI

from gradienteval import GradientEvaluator

def cleanup_and_exit(gaussian_process=None, error_status=0):
    """Cleanup function to ensure proper termination"""
    if gaussian_process is not None:
        try:
            # Get the process group id
            pgid = os.getpgid(gaussian_process.pid)
            
            # Force kill the entire process group immediately
            try:
                os.killpg(pgid, signal.SIGKILL)
            except OSError as e:
                print " * Error killing gaussian process group:", str(e)
            
        except OSError as e:
            print " * Error getting gaussian process group:", str(e)
    
    if MPI.Is_initialized() and not MPI.Is_finalized():
        try:
            MPI.Finalize()
        except:
            MPI.COMM_WORLD.Abort(error_status)
    
    sys.exit(error_status)

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()

# Define global variables
codeName = ""
dumpFileN = 'gau_EIn'  # gau-external.py dumps this file with data received from Gaussian
outPipeN = 'gau_EOut'  # Pipe for sending data back to Gaussian


if rank == 0:
    if len(sys.argv) < 2:
       print sys.argv
       print "Provide a .com file"
       sys.exit()
    
    gau_ComFile = sys.argv[1]
    inputPrefix = gau_ComFile.replace('.com','')
    print inputPrefix

    comFl = open(gau_ComFile, 'r') 

    inputTemplate = ''
    methodKey = ''
    userCs = []

    while True:
        line = comFl.readline()
        if line == '':
            break

        if line.find('!CODE ') > -1:
            L = line.upper().split()
            codeName = L[1]
            print 'Code name found:', codeName

        if ( line.upper().find('TEMPLATE') > -1 ) and ( line.upper().find(codeName) > -1 ) :
            line = comFl.readline()
            while line.upper().find('TEMPLATE') == -1:
                inputTemplate += line[1:] 
                #
                # MRCC-specific
                if line.upper().find('CALC=') > -1:
                    methodKey = line.rstrip().split('=')[1]
                #
                line = comFl.readline()

        if line.upper().find('SYMMETRY CONSTRAINTS') > -1:
            line = comFl.readline()
            while line.upper().find('SYMMETRY CONSTRAINTS') == -1:
                userCs.append(line)
                line = comFl.readline()

    comFl.close()

    if inputTemplate == '':
        print "No %s input template found" % codeName
        sys.exit()
    print " * Found %s input file template:" % codeName
    print inputTemplate

    if codeName == "":
        print "No codeName found" % codeName
        sys.exit()
 
    if codeName == "MRCC":
        if methodKey == '':
            print "MRCC method key is not found"
            sys.exit()
        print " * MRCC method: %s" % methodKey

    if userCs == []:
        print " ! No user-supplied symmetry constraints found; proceeding without."
    else:
        print " * User-supplied symmetry constraints found:"
        for line in userCs:
            print line,

gaussian_process = None
if rank == 0:
    # Launch gaussian.sh with better process management
    try:
        gaussian_process = subprocess.Popen(
            ["./gaussian.sh", gau_ComFile],
            shell=False,
            preexec_fn=os.setsid
        )
        print " * Rank %i %s launched gaussian.sh" % (rank, name)
    except OSError as e:  # Python 2.7 raises OSError instead of FileNotFoundError
        if e.errno == os.errno.ENOENT:
            print " * Error: gaussian.sh not found"
            cleanup_and_exit(error_status=127)
        elif e.errno == os.errno.EACCES:
            print " * Error: gaussian.sh not executable"
            cleanup_and_exit(error_status=126)
        else:
            print " * Error launching gaussian.sh: %s" % str(e)
            cleanup_and_exit(error_status=1)
    except Exception as e:
        print " * Error launching gaussian.sh: %s" % str(e)
        cleanup_and_exit(error_status=1)

# giving Gaussian time to do stuff
time.sleep(1.0)

comm.barrier()
# sending other ranks the external code name
codeName = comm.bcast(codeName, root=0)
# importing the apropriate module to be able to run estcode.executeCode()
try:
    exec 'import %s as estcode' % ( codeName.lower() ) 
except:
    moduleProblemMsg = "Cannot import module '%s'"
    raise Exception(moduleProblemMsg % codeName.lower() )

status = None 

while True:
    # Check gaussian.sh process status first
    if rank == 0 and gaussian_process is not None:
        status = gaussian_process.poll()
        if status is not None:
            if status == 0:
                print " * gaussian.sh completed successfully"
            else:
                print " * gaussian.sh failed with exit code %i" % status
                # Detailed error reporting based on exit code
                if status == 126:
                    print " * Error: gaussian.sh or g16 not executable"
                elif status == 127:
                    print " * Error: gaussian.sh or g16 not found in PATH"
                elif status == 137:
                    print " * Error: Process killed (possibly out of memory)"
                elif status == 1:
                    print " * Error: General gaussian.sh error (check output files)"
                elif status == -15:
                    print " * Error: Process terminated by SIGTERM"
                elif status == -9:
                    print " * Error: Process killed by SIGKILL"
                elif status < 0:
                    print " * Error: Process terminated by signal %i" % (-status)
                else:
                    print " * Error: Unknown error code %i" % status
            
            cleanup_and_exit(gaussian_process, error_status=abs(status))

    chunks = []

    #########################################################
    # Master process setting up stuff
    if rank == 0:
        # using a file instead of a pipe to get data from gau-external.py to avoid possible hang-ups
        print " * Rank %i %s waiting for Gaussian" % (rank, name)
        while not os.path.exists(dumpFileN):
            time.sleep(1.5)
            # Check gaussian.sh process status while waiting
            status = gaussian_process.poll()
            if status is not None:
                print " * Rank %i %s: gaussian.sh exited with status %i." % (rank, name, status)
                break 

    # Broadcast exit decision to all ranks
    comm.barrier()
    status = comm.bcast(status, root=0)

    # ALL ranks check exit status and cleanup if needed
    if status is not None:
        if rank == 0:
            print " * Rank 0 Initiating cleanup and exit on all ranks because Gaussinan quit with status %i " % status
        cleanup_and_exit(gaussian_process if rank == 0 else None, error_status=status)
                
    if rank == 0:
        thisMolecule = Molecule() 
        nAtoms, derivatives = thisMolecule.readGaussianExternal(dumpFileN) 
        os.remove(dumpFileN)
        print " * Rank %i %s coordinates from Gaussian received" % (rank, name)
        print thisMolecule.XYZonly()

        if userCs == []:
            print " * No symmetry constraints defined"
            thisGradientEv = GradientEvaluator(thisMolecule,filePrefix=inputPrefix,codeName=codeName)
        else:
            thisGradientEv = GradientEvaluator(thisMolecule,symmConstraints=userCs,filePrefix=inputPrefix,codeName=codeName)
            errMsg = thisGradientEv.parseSymmConstraints()
            if errMsg == "":
                print " * Parsed symmetry constraints:"
                print thisGradientEv.printSymmConstraints()
            else:
                print errMsg 
                cleanup_and_exit(gaussian_process, error_status=1)

        print " * Generating unique displacements"
        thisGradientEv.generateDisplacements() 
        fileList = thisGradientEv.generateESTinputs(inputTemplate)
        print " * Generated %i files" % len(fileList)

        for i in range(size): 
            chunks.append([])   

        i = 0
        for fileN in fileList: 
            chunks[i].append(fileN)
            i += 1
            if i == size: 
                i = 0 

    # Broadcast the chunks to all processes
    comm.barrier()
    fileList2Run = comm.scatter(chunks, root=0)

    error_info = []  # List to store error information

    for file2Run in fileList2Run:
        time.sleep(rank*0.05+0.01)
        
        print " > Rank %02i running %s for %s host: %s" % ( rank, codeName, file2Run, name )
        errcode = estcode.executeCode(file2Run,Verbose=False)
        if errcode != 0:
            #print " !!!\n !!! ERROR: Rank %02i  %s run for %s failed; host: %s\n !!!" % ( rank, codeName, file2Run, name )    
            # print " !!!\n !!! ERROR: Rank %02i %s run for %s failed; host: %s, status: %i\n !!!" % ( rank, codeName, file2Run, name, errcode )    
            error_info.append((rank, file2Run, name, errcode))
        else:
            print " > Rank %02i %s run for %s complete; host: %s, status: %i" % ( rank, codeName, file2Run, name, errcode )

    comm.barrier() # Make sure all runs are truly complete

    # Gather error information from all processes
    all_error_info = comm.gather(error_info, root=0)

    has_errors = False
    if rank == 0:
        # Master process checks for errors and prints error information
        print "\n * All %s runs finished.\n" % codeName

        for proc_error_info in all_error_info:
            if proc_error_info:
                has_errors = True
                for error_data in proc_error_info:
                    rank_num, file_name, host_name, error_code = error_data
                    # print " !!! Error on rank %02i: file: %s, host: %s, status: %i" % (rank_num, file_name, host_name, error_code)
                    print " !!! Error on rank %02i: file: %s, host: %s" % (rank_num, file_name, host_name)

    
    comm.barrier()
    # Broadcast exit decision to all ranks
    has_errors = comm.bcast(has_errors, root=0)
    
    comm.barrier()

    # ALL ranks check if there were errors and cleanup if needed
    if has_errors:
        if rank == 0:
            print " !!! Rank %i: Some %s runs failed, terminating..." % ( rank, codeName )
            print " !!! Rank %i: Initiating cleanup and exit on all ranks" % rank
            cleanup_and_exit(gaussian_process if rank == 0 else None, error_status=1)
        else:
            # print " * Rank %i %s received termination signal because some runs failed" % (rank, name)
            cleanup_and_exit(gaussian_process if rank == 0 else None, error_status=1)
    else:
        if rank == 0:
            print " * No errors reported in %s runs, proceeding... " % codeName

    if rank == 0:
        print " * Rank %i %s gathering data from %s runs" % (rank, name, codeName)

        zeroEnergyErrorMsg = thisGradientEv.readESToutputs(methodKey)

        if zeroEnergyErrorMsg != "":
            print zeroEnergyErrorMsg
            print " !!! Terminating..." 
            cleanup_and_exit(gaussian_process, error_status=1)

        print " * Rank %i %s evaluating gradients" % (rank, name)
        thisGradientEv.calculateGradients()
        print thisGradientEv.mol_.getAtomGradStr()

        print " * Rank %i %s sending data to Gaussian" % (rank, name)
        outPipe = open(outPipeN,"w")
        outPipe.write(thisGradientEv.mol_.GaussianExternalOutput())
        outPipe.close() 
        print " * Rank %i %s data to Gaussian sent" % (rank, name)

        # giving Gaussian time to work
        time.sleep(1.5)


    # Make sure all ranks are caught up
    comm.barrier()

    # Determine exit status (rank 0)
    if rank == 0:
        status = gaussian_process.poll()
        if status is not None:
            if status == 0:
                print " * gaussian.sh completed successfully"
            else:
                print " * gaussian.sh failed with exit code %i" % status
    else:
        status = None

    # Make sure all ranks are caught up
    comm.barrier()
    # Broadcast exit decision to all ranks
    status = comm.bcast(status, root=0)
    
    # ALL ranks check exit status and cleanup if needed
    if status is not None:
        if rank == 0:
            print " * Initiating cleanup and exit on all ranks"
        cleanup_and_exit(gaussian_process if rank == 0 else None, error_status=status)

    # If we're still here, return to main loop for next iteration












