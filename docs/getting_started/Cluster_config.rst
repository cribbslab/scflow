The .cgat.yml is placed in your home directory and when a pipeline is executed it will automatically prioritise the 
:file:`.cgat.yml` parameters over the cgatcore hard coded parameters. For example, adding the following to the
.cgat.yml file will implement cluster settings for PBSpro::


	memory_resource: mem
    
    options: -l walltime=00:10:00 -l select=1:ncpus=8:mem=1gb
    
    queue_manager: pbspro
    
    queue: NONE
    
    parallel_environment:  "dedicated"
    

 
