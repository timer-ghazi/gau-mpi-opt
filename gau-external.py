#!/usr/bin/env python
#
Header="""
# 
#
# --- Embarrassingly Parallel Gaussian Geometry Optimization with an External Code --- 
#                    High-Level Energies + Gaussian Optimizer 
#                        Q.T., Marquette U., 2018--2024
#                            
#                               gau-external.py 
#                          Gaussian Interface Script 
#                           
#
"""
# Is called by G16 via 'External' keyword
# * copies data from Gaussian to gau_EIn file for the main script to pick up
# * receives gradient/energy from the main script via gau_EOut pipe
# * dumps the data received into a file specified by Gaussian
# 
#
#       Thank you Markus for teaching me about pipes back in 2008!
#

import sys
import os
import os.path

def cleanup_files(files):
    """Safely remove files if they exist"""
    for f in files:
        try:
            if os.path.exists(f):
                os.remove(f)
        except OSError as e:
            print "Warning: Could not remove file %s: %s" % (f, str(e))

def write_file_safely(filename, content):
    """Safely write content to file"""
    try:
        with open(filename, "w") as f:
            f.write(content)
        return True
    except IOError as e:
        print "Error writing to %s: %s" % (filename, str(e))
        return False

def read_file_safely(filename):
    """Safely read from file"""
    try:
        with open(filename, "r") as f:
            return f.read()
    except IOError as e:
        print "Error reading from %s: %s" % (filename, str(e))
        return None

if __name__ == "__main__": 
    print Header
    DUMP_FILE = 'gau_EIn'
    PIPE_NAME = 'gau_EOut'
    
    # Validate command line arguments
    if len(sys.argv) < 4:
        print "Usage: %s <command> <input_file> <output_file>" % sys.argv[0]
        sys.exit(1)
        
    input_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Read Gaussian input
    gaussian_data = read_file_safely(input_file)
    if gaussian_data is None:
        sys.exit(1)
    
    # Write data for MPI script
    if not write_file_safely(DUMP_FILE, gaussian_data):
        cleanup_files([DUMP_FILE])
        sys.exit(1)
    
    # Create named pipe if needed
    try:
        if not os.path.exists(PIPE_NAME):
            os.mkfifo(PIPE_NAME)
    except OSError as e:
        print "Error creating pipe %s: %s" % (PIPE_NAME, str(e))
        cleanup_files([DUMP_FILE])
        sys.exit(1)
    
    # Read derivatives from pipe - simple blocking read
    print "gau-external.py: Waiting for derivatives..."
    try:
        with open(PIPE_NAME, "r") as pipe:
            derivatives = pipe.read()
    except IOError as e:
        print "Error reading from pipe %s: %s" % (PIPE_NAME, str(e))
        cleanup_files([DUMP_FILE])
        sys.exit(1)
    
    print "gau-external.py: Derivatives received:"
    print 
    
    # Write back to Gaussian
    if not write_file_safely(output_file, derivatives):
        cleanup_files([DUMP_FILE])
        print "gau-external.py: Derivatives sent"
        print 
        sys.exit(1)
    
    # Cleanup
    cleanup_files([DUMP_FILE])


