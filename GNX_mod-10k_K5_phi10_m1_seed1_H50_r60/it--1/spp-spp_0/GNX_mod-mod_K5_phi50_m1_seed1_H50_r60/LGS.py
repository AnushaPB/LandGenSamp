import geonomics as gnx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#set seed
np.random.seed(42)

#Make uniform array
def make_unif_array(n):
    """Makes a square array of ones, size n x n cells."""
    array = np.ones((n,n))
    return array
  
unifenv = make_unif_array(41)

#define default parameters (!these will be changed in the loop, I just define them here to create the variables!)
#Params to define
##K_factor = K
K = 0
##movement_distance_distr_param1 = md/dispersal_distance_distr_param1 = dd
m = 0
##phi = phi
phi = 0
##envs
env1 = unifenv
env2 = unifenv
  

                   #   ::::::          :::    :: :::::::::::#
             #::::::    ::::   :::      ::    :: :: ::::::::::: ::#
          #:::::::::     ::            ::   :::::::::::::::::::::::::#
        #::::::::::                      :::::::::: :::::: ::::::::  ::#
      #  : ::::  ::                    ::::  : ::    :::::::: : ::  :    #
     # GGGGG :EEEE: OOOOO   NN   NN   OOOOO   MM   MM IIIIII  CCCCC SSSSS #
    # GG     EE    OO   OO  NNN  NN  OO   OO  MM   MM   II   CC     SS     #
    # GG     EE   OO     OO NN N NN OO     OO MMM MMM   II   CC     SSSSSS #
    # GG GGG EEEE OO     OO NN  NNN OO     OO MM M MM   II   CC         SS #
    # GG   G EE    OO   OO  NN   NN  OO   OO  MM   MM   II   CC        SSS #
     # GGGGG :EEEE: OOOOO   NN   NN   OOOOO   MM   MM IIIIII  CCCCC SSSSS #
      #    : ::::::::               :::::::::: ::              ::  :   : #
        #:    :::::                    :::::: :::             :::::::  #
          #    :::                      :::::  ::              ::::: #
             #  ::                      ::::                      #
                   #                                        #
                      #  :: ::    :::             #


params = {
#--------------------------------------------------------------------------#

#-----------------#
#--- LANDSCAPE ---#
#-----------------#
    'landscape': {

    #------------#
    #--- main ---#
    #------------#
        'main': {
            #x,y (a.k.a. j,i) dimensions of the Landscape
            'dim':                      (41,41),
            #x,y resolution of the Landscape
            'res':                      (1,1),
            #x,y coords of upper-left corner of the Landscape
            'ulc':                      (0,0),
            #projection of the Landscape
            'prj':                      None,
            }, # <END> 'main'

    #--------------#
    #--- layers ---#
    #--------------#
        'layers': {

            #layer name (LAYER NAMES MUST BE UNIQUE!)
            'lyr_0': {

        #-------------------------------------#
        #--- layer num. 0: init parameters ---#
        #-------------------------------------#

                #initiating parameters for this layer
                'init': {

                    #parameters for a 'defined'-type Layer
                    'defined': {
                        #raster to use for the Layer
                        'rast':                   unifenv,
                        #point coordinates
                        'pts':                    None,
                        #point values
                        'vals':                   None,
                        #interpolation method {None, 'linear', 'cubic',
                        #'nearest'}
                        'interp_method':          None,

                        }, # <END> 'defined'

                    }, # <END> 'init'

                }, # <END> layer num. 0


            #layer name (LAYER NAMES MUST BE UNIQUE!)
            'lyr_1': {

        #-------------------------------------#
        #--- layer num. 1: init parameters ---#
        #-------------------------------------#

                #initiating parameters for this layer
                'init': {

                    #parameters for a 'defined'-type Layer
                    'defined': {
                        #raster to use for the Layer
                        'rast':                   env1,
                        #point coordinates
                        'pts':                    None,
                        #point values
                        'vals':                   None,
                        #interpolation method {None, 'linear', 'cubic',
                        #'nearest'}
                        'interp_method':          None,

                        }, # <END> 'defined'

                    }, # <END> 'init'

                }, # <END> layer num. 1
            
            #layer name (LAYER NAMES MUST BE UNIQUE!)
            'lyr_2': {

        #-------------------------------------#
        #--- layer num. 2: init parameters ---#
        #-------------------------------------#

                #initiating parameters for this layer
                'init': {

                    #parameters for a 'defined'-type Layer
                    'defined': {
                        #raster to use for the Layer
                        'rast':                   env2,
                        #point coordinates
                        'pts':                    None,
                        #point values
                        'vals':                   None,
                        #interpolation method {None, 'linear', 'cubic',
                        #'nearest'}
                        'interp_method':          None,

                        }, # <END> 'defined'

                    }, # <END> 'init'

                }, # <END> layer num. 2

    #### NOTE: Individual Layers' sections can be copy-and-pasted (and
    #### assigned distinct keys and names), to create additional Layers.


            } # <END> 'layers'

        }, # <END> 'landscape'


#-------------------------------------------------------------------------#
    
#-----------------#
#--- COMMUNITY ---#
#-----------------#
    'comm': {

        'species': {

            #species name (SPECIES NAMES MUST BE UNIQUE!)
            'spp_0': {

            #-----------------------------------#
            #--- spp num. 0: init parameters ---#
            #-----------------------------------#

                'init': {
                    #starting number of individs
                    'N':                250,
                    #carrying-capacity Layer name
                    'K_layer':          'lyr_0',
                    #multiplicative factor for carrying-capacity layer    
                    'K_factor':         K,                                
                    }, # <END> 'init'

            #-------------------------------------#
            #--- spp num. 0: mating parameters ---#
            #-------------------------------------#

                'mating'    : {
                    #age(s) at sexual maturity (if tuple, female first)
                    'repro_age':                0,
                    #whether to assign sexes
                    'sex':                      False,
                    #ratio of males to females
                    'sex_ratio':                1/1,
                    #whether P(birth) should be weighted by parental dist
                    'dist_weighted_birth':       False,
                    #intrinsic growth rate
                    'R':                        0.5,
                    #intrinsic birth rate (MUST BE 0<=b<=1)
                    'b':                        0.2,
                    #expectation of distr of n offspring per mating pair
                    'n_births_distr_lambda':    1,
                    #whether n births should be fixed at n_births_dist_lambda
                    'n_births_fixed':           True,
                    #ADDED BY AB: choose nearest mate
                    'choose_nearest_mate':      False,
                    #ADDED BY AB: choose nearest mate
                    'inverse_dist_mating':      False,
                    #radius of mate-search area         
                    'mating_radius':            1, #CHECK!
                    }, # <END> 'mating'

            #----------------------------------------#
            #--- spp num. 0: mortality parameters ---#
            #----------------------------------------#

                'mortality'     : {
                    #maximum age
                    'max_age':                      None,
                    #min P(death) (MUST BE 0<=d_min<=1)
                    'd_min':                        0,
                    #max P(death) (MUST BE 0<=d_max<=1)
                    'd_max':                        1,
                    #width of window used to estimate local pop density
                    'density_grid_window_width':    None,
                    }, # <END> 'mortality'

            #---------------------------------------#
            #--- spp num. 0: movement parameters ---#
            #---------------------------------------#

                'movement': {
                    #whether or not the species is mobile
                    'move':                                 True,
                    #mode of distr of movement direction
                    'direction_distr_mu':                   1,
                    #concentration of distr of movement direction
                    'direction_distr_kappa':                0,
                    #1st param of distr of movement distance  (mean)       
                    'movement_distance_distr_param1':       m, 
                    #2nd param of distr of movement distance (variance)
                    'movement_distance_distr_param2':       5e-8,
                    #movement distance distr to use ('levy' or 'wald')
                    'movement_distance_distr':              'levy',
                    #1st param of distr of dispersal distance
                    'dispersal_distance_distr_param1':      m,
                    #2nd param of distr of dispersal distance
                    'dispersal_distance_distr_param2':      5e-14,
                    #dispersal distance distr to use ('levy' or 'wald')
                    'dispersal_distance_distr':             'levy',
                    'move_surf'     : {
                        #move-surf Layer name
                        'layer':                'lyr_0',
                        #whether to use mixture distrs
                        'mixture':              True,
                        #concentration of distrs
                        'vm_distr_kappa':       12,
                        #length of approximation vectors for distrs
                        'approx_len':           5000,
                        }, # <END> 'move_surf'

                    },    # <END> 'movement'


            #---------------------------------------------------#
            #--- spp num. 0: genomic architecture parameters ---#
            #---------------------------------------------------#

                'gen_arch': {
                    #whether to jitter recomb bps, only needed to correctly track num_trees
                    'jitter_breakpoints':       False,
                    #file defining custom genomic arch
                    'gen_arch_file':            None,
                    #num of loci
                    'L':                        10000,
                    #num of chromosomes
                    'l_c':                      [100],
                    #starting allele frequency (None to draw freqs randomly)
                    'start_p_fixed':            0.5,
                    #whether to start neutral locus freqs at 0
                    'start_neut_zero':          False,
                    #genome-wide per-base neutral mut rate (0 to disable)
                    'mu_neut':                  0,
                    #genome-wide per-base deleterious mut rate (0 to disable)
                    'mu_delet':                 0,
                    #shape of distr of deleterious effect sizes
                    'delet_alpha_distr_shape':  0.2,
                    #scale of distr of deleterious effect sizes
                    'delet_alpha_distr_scale':  0.2,
                    #alpha of distr of recomb rates
                    'r_distr_alpha':            None,
                    #beta of distr of recomb rates
                    'r_distr_beta':             None,
                    #whether loci should be dominant (for allele '1')
                    'dom':                      False,
                    #whether to allow pleiotropy
                    'pleiotropy':               False,
                    #custom fn for drawing recomb rates
                    'recomb_rate_custom_fn':    None,
                    #number of recomb paths to hold in memory
                    'n_recomb_paths_mem':       int(1e4),
                    #total number of recomb paths to simulate
                    'n_recomb_paths_tot':       int(1e5),
                    #num of crossing-over events (i.e. recombs) to simulate
                    'n_recomb_sims':            10000,
                    #whether to generate recombination paths at each timestep
                    'allow_ad_hoc_recomb':      False,
                    #whether to save mutation logs
                    'mut_log':                  False,

                    'traits': {

                        #--------------------------#
                        #--- trait 1 parameters ---#
                        #--------------------------#
                        #trait name (TRAIT NAMES MUST BE UNIQUE!)
                        'trait_1': {
                            #trait-selection Layer name
                            'layer':                'lyr_1',
                            #polygenic selection coefficient    
                            'phi':                  phi,       
                            #number of loci underlying trait    
                            'n_loci':               4,          
                            #mutation rate at loci underlying trait
                            'mu':                   0,
                            #mean of distr of effect sizes
                            'alpha_distr_mu' :      0.25,
                            #variance of distr of effect size
                            'alpha_distr_sigma':    0,
                            #max allowed magnitude for an alpha value
                            'max_alpha_mag':        None,
                            #curvature of fitness function
                            'gamma':                1,
                            #whether the trait is universally advantageous
                            'univ_adv':             False
                            }, # <END> trait 0
                        
                        #--------------------------#
                        #--- trait 2 parameters ---#
                        #--------------------------#
                        #trait name (TRAIT NAMES MUST BE UNIQUE!)
                        'trait_2': {
                            #trait-selection Layer name
                            'layer':                'lyr_2',
                            #polygenic selection coefficient    
                            'phi':                  phi,       
                            #number of loci underlying trait    
                            'n_loci':               4,          
                            #mutation rate at loci underlying trait
                            'mu':                   0,
                            #mean of distr of effect sizes
                            'alpha_distr_mu' :      0.25,
                            #variance of distr of effect size
                            'alpha_distr_sigma':    0,
                            #max allowed magnitude for an alpha value
                            'max_alpha_mag':        None,
                            #curvature of fitness function
                            'gamma':                1,
                            #whether the trait is universally advantageous
                            'univ_adv':             False
                            }, # <END> trait 0


    #### NOTE: Individual Traits' sections can be copy-and-pasted (and
    #### assigned distinct keys and names), to create additional Traits.


                        }, # <END> 'traits'

                    }, # <END> 'gen_arch'


                }, # <END> spp num. 0



    #### NOTE: individual Species' sections can be copy-and-pasted (and
    #### assigned distinct keys and names), to create additional Species.


            }, # <END> 'species'

        }, # <END> 'comm'


#------------------------------------------------------------------------#

#-------------#
#--- MODEL ---#
#-------------#
    'model': {
        #total Model runtime (in timesteps)
        'T':            1001,
        #min burn-in runtime (in timesteps)
        'burn_T':       100,
        #seed number
        'num':          42,
        #time step interval for simplication of tskit tables
        'tskit_simp_interval':      100,


        #-----------------------------#
        #--- iterations parameters ---#
        #-----------------------------#
        'its': {
            #num iterations
            'n_its':            1,
            #whether to randomize Landscape each iteration
            'rand_landscape':   False,
            #whether to randomize Community each iteration
            'rand_comm':        False,
            #whether to burn in each iteration
            'repeat_burn':      False,
            }, # <END> 'iterations'
      
        #-----------------------------------#
        #--- data-collection parameters ---#
        #-----------------------------------#
        'data': {
            'sampling': {
                #sampling scheme {'all', 'random', 'point', 'transect'}
                'scheme':               'all',
                #when to collect data
                'when':                 100,
                #whether to save current Layers when data is collected
                'include_landscape':    False,
                #whether to include fixed loci in VCF files
                'include_fixed_sites':  True,
                },
            'format': {
                #format for genetic data {'vcf', 'fasta'}
                'gen_format':           ['vcf', 'fasta'],
                #format for vector geodata {'csv', 'shapefile', 'geojson'}
                'geo_vect_format':      'csv',
                #format for raster geodata {'geotiff', 'txt'}
                'geo_rast_format':      'geotiff',
                },
            }, #<END> 'data'
      
        #-----------------------------------#
        #--- stats-collection parameters ---#
        #-----------------------------------#
        'stats': {
            #number of individs at time t
            'Nt': {
                #whether to calculate
                'calc':     True,
                #calculation frequency (in timesteps)
                'freq':     1,
                },
            #heterozgosity
            'het': {
                #whether to calculate
                'calc':     True,
                #calculation frequency (in timesteps)
                'freq':     5,
                #whether to mean across sampled individs
                'mean': False,
                },
            #minor allele freq
            'maf': {
                #whether to calculate
                'calc':     True,
                #calculation frequency (in timesteps)
                'freq':     5,
                },
            #mean fitness
            'mean_fit': {
                #whether to calculate
                'calc':     True,
                #calculation frequency (in timesteps)
                'freq':     5,
                },
            #linkage disequilibirum
            'ld': {
                #whether to calculate
                'calc':     False,
                #calculation frequency (in timesteps)
                'freq':     100,
                },
            }, # <END> 'stats'

        } # <END> 'model'

    } # <END> params


#TESTING LOOP WILL PARALLELIZE LATER
for K in [5]:
  for phi in [0.5]:
    for m in[1]:
      for seed in [1]:
        for H in [0.5]:
          for r in [0.6]: 
            #directory where input/output data will be stored (scratch)
            dir = "/global/scratch/anushabishop/LGS/"
            #note: currently gnx dumps output files where the script is run
            
            #get env layers
            env1 = np.genfromtxt(dir+"MNLM/seed"+str(seed)+"_env1_H"+str(int(H*100))+"_r"+str(int(r*100))+".csv", delimiter=',')
            env2 = np.genfromtxt(dir+"MNLM/seed"+str(seed)+"_env2_H"+str(int(H*100))+"_r"+str(int(r*100))+".csv", delimiter=',')
            
            #redefine params
            params['landscape']['layers']['lyr_1']['init']['defined']['rast'] = env1
            params['landscape']['layers']['lyr_2']['init']['defined']['rast'] = env2
            params['comm']['species']['spp_0']['init']['K_factor'] = K
            params['comm']['species']['spp_0']['movement']['movement_distance_distr_param1'] = m
            params['comm']['species']['spp_0']['movement']['dispersal_distance_distr_param1'] = m
            params['comm']['species']['spp_0']['gen_arch']['traits']['trait_1']['phi'] = phi
            
            #print params to confirm proper params were used (in output)
            print(params)
            
            #make our params dict into a proper Geonomics ParamsDict object
            mod_name = "mod_K"+str(K)+"_phi"+str(int(phi*100))+"_m"+str(m)+"_seed"+str(seed)+"_H"+str(int(H*100))+"_r"+str(int(r*100))
            print(mod_name)
            params = gnx.make_params_dict(params, mod_name)
            #then use it to make a model
            mod = gnx.make_model(parameters=params, verbose=True)
            
            #burn-in model
            mod.walk(T=10000, mode='burn')
                                                                                                                           
            #save and print all of the non-neutral loci
            loci_df = pd.DataFrame()
            loci_df['trait1'] = mod.comm[0].gen_arch.traits[0].loci
            loci_df['trait2'] = mod.comm[0].gen_arch.traits[1].loci
            loci_df.to_csv(dir+"nnloci_"+mod_name+".csv")  
            print("\nNON-NEUTRAL LOCI:")
            print(mod.comm[0].gen_arch.nonneut_loci)
                                                                                                                           
            #run the model for 1000 steps
            mod.walk(1000)
            