# MPI-Multithread-Framework
3 sequential algorithm (simplex, game-of-life and heatmap) were sped up using using MPI Framework

Full specification of those algorithms and constraints about their speeding up can be seen [here](https://github.com/mdodovic/MPI-Multithread-Framework/blob/main/problems_description.pdf).

## [task1](https://github.com/mdodovic/MPI-Multithread-Framework/tree/main/task1_simplex)

Sequential code is modified and then both basic and revised codes are sped up by dividing all of the _unit_of_works_ among the processes.

## [task2](https://github.com/mdodovic/MPI-Multithread-Framework/tree/main/task2_gameoflife)

Sequential code is sped up by dividing the complete work into processes. Some of the data (the first and the last row of the processes' part of the matrix) are exchanged between processes with adjacent rank. 

## [task3](https://github.com/mdodovic/MPI-Multithread-Framework/tree/main/task3_hotspot/data)

Sequential code is sped up by dividing the complete work into processes. Some of the data (the first and the last row of the processes' part of the matrix) are exchanged between processes with adjacent rank. 

## [task4](https://github.com/mdodovic/MPI-Multithread-Framework/tree/main/task4_hotspot/data)

Sequential code is sped up using the _Manager-Worker_ technique. One (rank = 0) process was declared as master and all of the others processes requested the job to be done. After the job was done, the result was sent to the master.  
