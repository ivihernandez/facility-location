#standard imports
import time
import gc
import os
#non standard imports

#ivan's imports
import runner

        
def main():
    
    
    experimentsFolder = r'.\tests-folder'
    for element in os.listdir(experimentsFolder):
        
        path = os.path.join(experimentsFolder,element)
        if os.path.isdir(path):
            #the element is a folder, search it and try to run them
            files = os.listdir(path)
            for myfile in files:
                myfilepath = os.path.join(path, myfile)
                print "the parameter will be ", myfilepath
                runner.run_experiment(myfilepath)
                gc.collect()
        elif os.path.isfile(path):
            print "the parameter will be ",path
            runner.run_experiment(path)
            gc.collect()
        else:
            print "is nothing"
   
if __name__ == '__main__':
    print 'Program started'
    start = time.clock()
    main()
    end = time.clock()
    print 'Program finished in ',(end-start)/60.0 ,' minutes'
    