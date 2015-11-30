#include "NESTNodeSynapse.h"
#include "H5SynapseLoader.h"

#include "timer/stopwatch.h"


#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <fstream>
#include <sstream>

#include <omp.h>
#include <cmath>

#include <limits>


int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  
  int RANK, NUM_PROCESSES;
  MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
  MPI_Comm_size(MPI_COMM_WORLD, &NUM_PROCESSES);
  
  
  std::string syn_filename = "/gpfs/bbp.cscs.ch/scratch/gss/bgq/schumann/syn_long_full_20151112.h5";
  
  
  uint64_t n_readSynapses =0;
  uint64_t n_SynapsesInDatasets =0;
  uint64_t n_memSynapses =0;
  
  int c;
  
  while ((c = getopt (argc, argv, "f:")) != -1)
    switch (c)
    {
      case 'f':
	  syn_filename = optarg;
	break;
      default:
	break;
    }
  
  H5SynapsesLoader synloader(syn_filename,n_readSynapses,n_SynapsesInDatasets);
  
  
  
  uint64_t b = -1;
  uint64_t e = -1;
  
  
  while ((c = getopt (argc, argv, "b:e:")) != -1)
    switch (c)
    {
      case 'b':
	  b = (uint64_t)atoi(optarg)*(uint64_t)1e6*(uint64_t)1e4;
	  synloader.setFirstSyn(b);
	break;
      case 'e':
	  e = (uint64_t)atoi(optarg)*(uint64_t)1e6*(uint64_t)1e4;
	  synloader.setLastSyn(e);
	break;
      case 'f':
	  syn_filename = optarg;
	break;
      default:
	break;
    }
  
  
  //number of synapses per iteration effects memory consumption and speed of the import module
  uint64_t nos = 1e7; 
  
  
  //std::vector < uint64_t > num_pre_synapses(74107932,0);
  
  //std::vector < uint64_t > num_pre_synapses(74107932,0);
  
 
  double avg_delay=0;
  double avg_weight=0;
  double avg_U0=0;
  double avg_TauRec=0;
  double avg_TauFac=0;
  
  float max_delay=std::numeric_limits<float>::min();
  float max_weight=std::numeric_limits<float>::min();
  float max_U0=std::numeric_limits<float>::min();
  float max_TauRec=std::numeric_limits<float>::min();
  float max_TauFac=std::numeric_limits<float>::min();
  
  float min_delay=std::numeric_limits<float>::max();
  float min_weight=std::numeric_limits<float>::max();
  float min_U0=std::numeric_limits<float>::max();
  float min_TauRec=std::numeric_limits<float>::max();
  float min_TauFac=std::numeric_limits<float>::max();
  
  
    
    
  
  uint64_t it =0;
  
  
  NESTSynapseList synapses_;
  
  
  
  //load datasets from files
  while (!synloader.eof())
  {
    nest::Stopwatch::timestamp_t before_io=nest::Stopwatch::get_timestamp();
    synloader.iterateOverSynapsesFromFiles(synapses_, nos);   
    
    nest::Stopwatch::timestamp_t dur_io=nest::Stopwatch::get_timestamp()-before_io;

    std::cout << "run\t" << RANK << "\t" << it << "\t" << synapses_.size() << "\t" << n_memSynapses << "\t" << dur_io /1000 << std::endl;
    
    n_memSynapses+=synapses_.size();
    
    //#pragma omp parallel 
    //{
    //  const int thread = omp_get_thread_num();
    //  const int num_threads = omp_get_num_threads();
      
    //  const double fact = (double)74107932/(double)num_threads;
      
    //  const int start = round(thread*fact);
    //  const int end = round((thread+1)*fact); 
      
      for (int i=0; i<synapses_.size(); i++) {
	//const int tn = synapses_[i].target_neuron_;
	//num_pre_synapses[tn]++;
	avg_delay += (synapses_[i].delay - avg_delay) / (i+1);
	avg_weight += (synapses_[i].weight - avg_weight) / (i+1);
	avg_U0 += (synapses_[i].U0 - avg_U0) / (i+1);
	avg_TauRec += (synapses_[i].TauRec - avg_TauRec) / (i+1);
	avg_TauFac += (synapses_[i].TauFac - avg_TauFac) / (i+1);
	
	max_delay = std::max(max_delay, synapses_[i].delay);
	max_weight = std::max(max_weight, synapses_[i].weight);
	max_U0 = std::max(max_U0, synapses_[i].U0);
	max_TauRec = std::max(max_TauRec, synapses_[i].TauRec);
	max_TauFac = std::max(max_TauFac, synapses_[i].TauFac);
	
	min_delay = std::min(min_delay, synapses_[i].delay);
	min_weight = std::min(min_weight, synapses_[i].weight);
	min_U0 = std::min(min_U0, synapses_[i].U0);
	min_TauRec = std::min(min_TauRec, synapses_[i].TauRec);
	min_TauFac = std::min(min_TauFac, synapses_[i].TauFac);
	
      }
    //}
    
    //synapses_.clear();
    
    it++;
  }
  
  std::cout << "DEBUG\tEND\t" << RANK << std::endl; 
  
  
  //std::vector < uint64_t > global_num_pre_synapses(74107932,0);
  //MPI_Allreduce(&num_pre_synapses[0], &global_num_pre_synapses[0], num_pre_synapses.size(), MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  
  double global_avg_delay;
  double global_avg_weight;
  double global_avg_U0;
  double global_avg_TauRec;
  double global_avg_TauFac;
  
  float global_max_delay;
  float global_max_weight;
  float global_max_U0;
  float global_max_TauRec;
  float global_max_TauFac;
  
  float global_min_delay;
  float global_min_weight;
  float global_min_U0;
  float global_min_TauRec;
  float global_min_TauFac;
 
  MPI_Allreduce(&avg_delay, &global_avg_delay, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&avg_weight, &global_avg_weight, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&avg_U0, &global_avg_U0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&avg_TauRec, &global_avg_TauRec, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&avg_TauFac, &global_avg_TauFac, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  global_avg_delay/=NUM_PROCESSES;
  global_avg_weight/=NUM_PROCESSES;
  global_avg_U0/=NUM_PROCESSES;
  global_avg_TauRec/=NUM_PROCESSES;
  global_avg_TauFac/=NUM_PROCESSES;
  
  MPI_Allreduce(&max_delay, &global_max_delay, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&max_weight, &global_max_weight, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&max_U0, &global_max_U0, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&max_TauRec, &global_max_TauRec, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&max_TauFac, &global_max_TauFac, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  
  MPI_Allreduce(&min_delay, &global_min_delay, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&min_weight, &global_min_weight, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&min_U0, &global_min_U0, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&min_TauRec, &global_min_TauRec, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&min_TauFac, &global_min_TauFac, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  
  
  if (RANK==0){
    /*for (int i=0; i<global_num_pre_synapses.size(); i++) {
      std::cout << "output\t" <<  i << "\t" <<  global_num_pre_synapses[i] << "\n";
    }*/
    /*std::stringstream ss;
    ss << "output_" << b << "_" << e << ".bin";
    std::ofstream myFile (ss.str().c_str(), std::ios::out | std::ios::binary);
    
    myFile.write ((char*)&global_num_pre_synapses[0], global_num_pre_synapses.size()*sizeof(uint64_t));
    
    myFile.close();*/
    
    std::cout << "delay \taverage=" << global_avg_delay << "\tmax=" << global_max_delay << "\tmin=" << global_min_delay << std::endl;
    std::cout << "weight \taverage=" << global_avg_weight << "\tmax=" << global_max_weight << "\tmin=" << global_min_weight << std::endl;
    std::cout << "U0 \taverage=" << global_avg_U0 << "\tmax=" << global_max_U0 << "\tmin=" << global_min_U0 << std::endl;
    std::cout << "TauRec \taverage=" << global_avg_TauRec << "\tmax=" << global_max_TauRec << "\tmin=" << global_min_TauRec << std::endl;
    std::cout << "TauFac \taverage=" << global_avg_TauFac << "\tmax=" << global_max_TauFac << "\tmin=" << global_min_TauFac << std::endl;
    
    
  }
  
  std::cout << "validation\t" << RANK << "\t" << n_readSynapses << "\t" << n_SynapsesInDatasets << "\t" << n_memSynapses << "\n";


  MPI_Finalize();
  return 0;

}