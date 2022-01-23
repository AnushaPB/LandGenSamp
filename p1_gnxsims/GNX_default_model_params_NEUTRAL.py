# GNX_default_model_params_NEUTRAL.py

# This is a parameters file generated by Geonomics
# (by the gnx.make_parameters_file() function).


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
#-----------------------------------------------------------------------------#

#-----------------#
#--- LANDSCAPE ---#
#-----------------#
    'landscape': {

    #------------#
    #--- main ---#
    #------------#
        'main': {
            #x,y (a.k.a. j,i) dimensions of the Landscape
            'dim':                      (20,20),
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

                    #parameters for a 'random'-type Layer
                    'random': {
                        #number of random points
                        'n_pts':                        500,
                        #interpolation method {'linear', 'cubic', 'nearest'}
                        'interp_method':                'linear',

                        }, # <END> 'random'

                    }, # <END> 'init'

                }, # <END> layer num. 0



    #### NOTE: Individual Layers' sections can be copy-and-pasted (and
    #### assigned distinct keys and names), to create additional Layers.


            } # <END> 'layers'

        }, # <END> 'landscape'


#-----------------------------------------------------------------------------#

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
                    'K_factor':         1,
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
                    #intrinsic growth rate
                    'R':                        0.5,
                    #intrinsic birth rate (MUST BE 0<=b<=1)
                    'b':                        0.2,
                    #expectation of distr of n offspring per mating pair
                    'n_births_distr_lambda':    1,
                    #whether n births should be fixed at n_births_dist_lambda
                    'n_births_fixed':           True,
                    #radius of mate-search area (None, for panmixia)
                    'mating_radius':            10,
                    #whether individs should choose nearest neighs as mates
                    'choose_nearest_mate':        False,
                    #whether mate-choice should be inverse distance-weighted
                    'inverse_dist_mating':      False,
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
                    'direction_distr_mu':                   0,
                    #concentration of distr of movement direction
                    'direction_distr_kappa':                0,
                    #1st param of distr of movement distance
                    'movement_distance_distr_param1':       0.01,
                    #2nd param of distr of movement distance
                    'movement_distance_distr_param2':       0.5,
                    #movement distance distr to use ('lognormal','levy','wald')
                    'movement_distance_distr':              'lognormal',
                    #1st param of distr of dispersal distance
                    'dispersal_distance_distr_param1':      -1,
                    #2nd param of distr of dispersal distance
                    'dispersal_distance_distr_param2':      0.05,
                    #dispersal distance distr to use ('lognormal','levy','wald')
                    'dispersal_distance_distr':             'lognormal',
                    },    # <END> 'movement'


            #---------------------------------------------------#
            #--- spp num. 0: genomic architecture parameters ---#
            #---------------------------------------------------#

                'gen_arch': {
                    #file defining custom genomic arch
                    'gen_arch_file':            None,
                    #num of loci
                    'L':                        100,
                    #fixed starting allele freq; None/False -> rand; True -> 0.5
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
                    'r_distr_alpha':            0.5,
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
                    'n_recomb_sims':            10_000,
                    #whether to generate recombination paths at each timestep
                    'allow_ad_hoc_recomb':      False,
                    #whether to jitter recomb bps, to correctly track num_trees
                    'jitter_breakpoints':       False,
                    #whether to save mutation logs
                    'mut_log':                  False,
                    #whether to use tskit (to record full spatial pedigree)
                    'use_tskit':                True,
                    #time step interval for simplication of tskit tables
                    'tskit_simp_interval':      100,


                    }, # <END> 'gen_arch'


                }, # <END> spp num. 0



    #### NOTE: individual Species' sections can be copy-and-pasted (and
    #### assigned distinct keys and names), to create additional Species.


            }, # <END> 'species'

        }, # <END> 'comm'


#-----------------------------------------------------------------------------#

#-------------#
#--- MODEL ---#
#-------------#
    'model': {
        #total Model runtime (in timesteps)
        'T':            100,
        #min burn-in runtime (in timesteps)
        'burn_T':       30,
        #seed number
        'num':          None,


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
            #whether to randomize GenomicArchitectures each iteration
            'rand_genarch':     True,
            #whether to burn in each iteration
            'repeat_burn':      False,
            }, # <END> 'iterations'



        } # <END> 'model'

    } # <END> params
