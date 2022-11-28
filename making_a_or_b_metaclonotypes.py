#!/usr/bin/env python3

"""
This script was adapted from the TCRdist3 documentation. 
I revised it to work for our data and to do our filtering.
I also wrote a script to sample a reference alpha chain repertoire to compare against. 
"""
import sys
import os
import numpy as np
import pandas as pd
import scipy.sparse
import matplotlib.pyplot as plt
from tcrdist.paths import path_to_base
from tcrdist.repertoire import TCRrep
from tcrdist.automate import auto_pgen
from tcrsampler.sampler import TCRsampler
from tcrdist.background import get_stratified_gene_usage_frequency
from tcrdist.background import sample_britanova
from tcrdist.background import make_gene_usage_counter, get_gene_frequencies
from tcrdist.background import calculate_adjustment, make_gene_usage_counter
from tcrdist.background import make_vj_matched_background, make_flat_vj_background
from tcrdist.ecdf import distance_ecdf, _plot_manuscript_ecdfs
from tcrdist.centers import calc_radii
from tcrdist.public import _neighbors_variable_radius
from tcrdist.public import _neighbors_sparse_variable_radius
from tcrdist.regex import _index_to_regex_str
from general_sample_britanova_fnct import sample_df_bkgd
    ## This is a script that I wrote.

def find_metaclonotypes(
    user_chain,
    input_file,
    out_path,
    ncpus = 4, 
    seed = 3434,
    min_subjects= 2,
    min_redundancy = 1):
    """
    This functions encapsulates a complete 
    workflow for finding meta-clonotypes in antigen-enriched data.
    min_subjects : int 
        Minimum publicity of Metaclonotype
    min_redundancy : int
        minimum non-redundant TCRs in lower level metaclonotypes
    """
    np.random.seed(seed)
    ############################################################################
    # Step 1: Select and load a antigen-enriched (sub)repertoire.           ####
    ############################################################################

    antigen_enriched_name= "stim_and_ref_input.csv"
    print(f"INITIATING A TCRrep() with", antigen_enriched_name)

    if(user_chain == "beta"):
        input_df = pd.read_csv(input_file, low_memory= False).astype(str)
            ## Doing this b/c pandas was treating my .csv weirdly 
        input_df = input_df[(input_df['v_b_gene'].notna()) &
                            (input_df['j_b_gene'].notna()) &
                            (~input_df.cdr3_b_aa.str.contains('n'))  & ## Invalid aa seqs removed
                            (~input_df.cdr3_b_aa.str.contains('\*')) ]
        input_df = input_df.sample(replace= False, n= 10001)
            ## Subsampling returns pw_beta, but not subsampling doesn't?
            ## This has to do w/ how TCRrep works it gives a warning when it should give an error. It doesn't compute distances for large datasets.
        # input_df= input_df[(input_df['v_b_gene'].notna()) & 
        #                    (input_df['j_b_gene'].notna()) &
        #                    (input_df['cdr3_b_aa'] != "" ) &
        #                    (input_df['v_b_gene']  != "" ) &
        #                    (input_df['j_b_gene']  != "" )] ## Are empty entries an issue?
        input_df= input_df[['subject','antigen','v_b_gene','cdr3_b_aa',"j_b_gene"]]

    elif(user_chain == "alpha"):
        input_df = pd.read_csv(input_file, low_memory= False).astype(str)
        input_df = input_df[(input_df['v_a_gene'].notna()) &
                            (input_df['j_a_gene'].notna()) &
                            (~input_df.cdr3_a_aa.str.contains('n'))  & ## Invalid aa seqs removed
                            (~input_df.cdr3_a_aa.str.contains('\*')) &
                            (~input_df.cdr3_a_aa.str.contains("l")) &
                            (~input_df.cdr3_a_aa.str.contains("#")) &
                            (~input_df.cdr3_a_aa.str.contains("U")) &
                            (~input_df.cdr3_a_aa.str.contains("O")) &
                            (~input_df.cdr3_a_aa.str.contains("f")) &
                            (~input_df.cdr3_a_aa.str.contains("k")) &
                            (~input_df.cdr3_a_aa.str.contains("\(")) &
                            (~input_df.cdr3_a_aa.str.contains("1")) &
                            (~input_df.cdr3_a_aa.str.contains("4")) &
                            (~input_df.cdr3_a_aa.str.contains("c")) &
                            (~input_df.cdr3_a_aa.str.contains("o")) &
                            (~input_df.cdr3_a_aa.str.contains("J"))]

        input_df= input_df[['subject','antigen','v_a_gene','cdr3_a_aa','j_a_gene']]

    from tcrdist.repertoire import TCRrep
    tr = TCRrep(cell_df = input_df, 
                organism = "human", 
                chains = [user_chain], 
                compute_distances = True, ## True produces pw_beta only if less than 10000 rows
                cpus= ncpus)
               
    if input_df.shape[0] > 10000: 
        if user_chain == "beta":
            print("CALCULATING pw_beta for large datasets")
            tr.compute_sparse_rect_distances(radius=50, chunk_size=100)
            tr.pw_beta = tr.rw_beta.copy()
            #tr.pw_beta = tr.pw_beta.toarray()
            #print(tr.pw_beta.nbytes)
            
            # print(tr.pw_beta[1:10])
            # print(dir(tr))
            # quit()
        elif user_chain == "alpha":
            print("CALCULATING pw_alpha for large datasets")
            tr.compute_sparse_rect_distances(radius=50, chunk_size=100)
            tr.pw_alpha = tr.rw_alpha.copy()
            #tr.pw_alpha = tr.pw_alpha.toarray()
            #print(tr.pw_alpha.nbytes) 
            
    # print(tr.show_incomplete().columns)
    # print(tr.show_incomplete()['v_b_gene'].unique())
    # tr.show_incomplete().to_csv("/stor/work/Ehrlich_COVID19/Users/John/projects/stim_TCRs/TCRdist/tcrdist3_analysis/output/incomplete.csv")
    tr.cpus = ncpus
    tr.clone_df['radius'] = 50

    

    from tcrdist.automate import auto_pgen
    print(f"COMPUTING PGEN WITH OLGA (Sethna et al 2018)")
    print("FOR ANTIGEN-ENRICHED CLONES TO BE USED FOR SUBSEQUENT ANALYSES")
    auto_pgen(tr)

    ###############################################################################
    #### Splitting function now b/c tcrdist code base has the chains hardcoded ####
    ###############################################################################

    if user_chain == "beta":
        print(f"USING tcrsampler TO CONSTRUCT A CUSTOM V-J MATCHED BACKGROUND")
        from tcrsampler.sampler import TCRsampler
        ts = TCRsampler(default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv')
        # Stratify sample so that each subject contributes similarly to estimate of 
        # gene usage frequency
        from tcrdist.background import get_stratified_gene_usage_frequency
        ts = get_stratified_gene_usage_frequency(ts = ts, replace = True) 
        # Synthesize an inverse probability weighted V,J gene background that matches 
        # usage in your enriched repertoire 
        df_vj_background = tr.synthesize_vj_matched_background(ts = ts, chain = user_chain)
        # Get a randomly drawn stratified sampler of beta, cord blood from 
        # Britanova et al. 2016 
        # Dynamics of Individual T Cell Repertoires: From Cord Blood to Centenarians
        from tcrdist.background import  sample_britanova

        df_britanova_100K = sample_britanova(size = 100000)
        # Append frequency columns using, using sampler above
        df_britanova_100K = get_gene_frequencies(ts = ts, 
                                                 df = df_britanova_100K)
        df_britanova_100K['weights'] = 1
        df_britanova_100K['source']  = "stratified_random"
    
        # Combine the two parts of the background into a single DataFrame
        df_bkgd = pd.concat([df_vj_background.copy(),
                             df_britanova_100K.copy()], axis = 0).\
                             reset_index(drop = True)  
        background_outfile = os.path.join(out_path, 
                                          f"{antigen_enriched_name}_{user_chain}.olga100K_brit100K_bkgd.csv")
        print(f'WRITING {background_outfile}')
        df_bkgd.to_csv(background_outfile, index = False)
        # Load the background to a TCRrep without computing pairwise distances 
        # (i.e., compute_distances = False)
        tr_bkgd = TCRrep(
            cell_df = df_bkgd,
            organism = "human", 
            chains = [user_chain], 
            compute_distances = False)
    
        ############################################################################
        # Step 4: Calculate Distances                                          #####
        ############################################################################
        print(f"COMPUTING RECTANGULARE DISTANCE")
        # print('pw_beta' in dir(tr))
        # print('rw_beta' in dir(tr))
        # print('pw_beta' in dir(tr) and 'rw_beta' in dir(tr))
        # print("Entering function")
        tr.compute_sparse_rect_distances(
            df = tr.clone_df,  ## commenting out b/c changing compute_distances to false, this will produce rw_beta output
            df2 = tr_bkgd.clone_df,
            radius=50,
            chunk_size = 100)
            ## This produces rw_beta

        # print('pw_beta' in dir(tr))
        # print('rw_beta' in dir(tr))
        # print('pw_beta' in dir(tr) and 'rw_beta' in dir(tr))

        # print(tr.rw_beta.get_shape())
        # print(tr.pw_beta.size)
            ## They don't have the same shape. 

        scipy_path= os.path.join(out_path, 
                                f"{antigen_enriched_name}_{user_chain}.rw_beta.npz") 
        scipy.sparse.save_npz(scipy_path, tr.rw_beta)
            # Tip: For larger dataset you can use a sparse implementation: 
            # 30.8 s ± 0 ns per loop ; tr.cpus = 6
            # %timeit -r tr.compute_sparse_rect_distances(df = tr.clone_df, df2 = tr_bkdg.clone_df,radius=50, chunk_size=85)
        ############################################################################
        # Step 5: Examine Density ECDFS                                        #####
        ############################################################################
            # Investigate the density of neighbors to each TCR, based on expanding 
            # distance radius.
        from tcrdist.ecdf import distance_ecdf, _plot_manuscript_ecdfs
        import matplotlib.pyplot as plt
            # Compute empirical cumulative density function (ecdf)
            # Compare Antigen Enriched TCRs (against itself).
        #print(dir(tr))
        #print(type(tr.pw_beta))
        thresholds, antigen_enriched_ecdf = distance_ecdf(
            tr.pw_beta, ## By the default this is tr.pw_beta, but when I make the input df bigger it doesn't output .pw_beta, even if I set compute_distances to True. This feels like a github issue that I should raise at some point. I'll need to produce a reproducilbe example. 
            thresholds=range(0,50,2))
            # Compute empirical cumulative density function (ecdf)
            # Compare Antigen Enriched TCRs (against) 200K probability 
            # inverse weighted background
        thresholds, background_ecdf = distance_ecdf(
            tr.rw_beta,
            thresholds=range(0,50,2),
            weights= tr_bkgd.clone_df['weights'], 
            absolute_weight = True)
            # plot_ecdf similar to tcrdist3 manuscript #
        antigen_enriched_ecdf[antigen_enriched_ecdf == antigen_enriched_ecdf.min()] = 1E-10
        f1 = _plot_manuscript_ecdfs(
            thresholds, 
            antigen_enriched_ecdf, 
            ylab= 'Proportion of Antigen Enriched TCRs', 
            cdr3_len=tr.clone_df.cdr3_b_aa.str.len(), 
            min_freq=1E-10)
        f1.savefig(os.path.join(out_path, f'{antigen_enriched_name}_{user_chain}.ecdf_AER_plot.png'))
        f2 = _plot_manuscript_ecdfs(
            thresholds,
            background_ecdf,
            ylab= 'Proportion of Reference TCRs',
            cdr3_len=tr.clone_df.cdr3_b_aa.str.len(),
            min_freq=1E-10)
        f2.savefig(os.path.join(out_path, f'{antigen_enriched_name}_{user_chain}.ecdf_BUR_plot.png'))

        ############################################################################
        # Step 6: Find optimal radii  (theta = 1E5                             #####
        ############################################################################
        # To ascertain which meta-clonotypes are likely to be most specific, 
        # take advantage of an existing function <bkgd_cntrl_nn2>.                                                                                                                                  

        level_tag = '1E5'
        from tcrdist.neighbors import bkgd_cntl_nn2
        centers_df  = bkgd_cntl_nn2(
            tr               = tr,
            tr_background    = tr_bkgd,
            weights          = tr_bkgd.clone_df.weights,
            ctrl_bkgd        = 10**-5, 
            col              = 'cdr3_b_aa',
            add_cols         = ['v_b_gene', 'j_b_gene'],
            ncpus            = ncpus,
            include_seq_info = True,
            thresholds       = [x for x in range(0,50,2)],
            generate_regex   = True,
            test_regex       = True,
            forced_max_radius = 36)
            ## I'm pretty sure this is the function that caused the MemoryError

        ############################################################################
        # Step 6.2: (theta = 1E5) ALL meta-clonotypes .tsv file                   ##
        ############################################################################
        # save center to project_path for future use
        centers_df.to_csv(os.path.join(out_path, 
                                       f'{antigen_enriched_name}_{user_chain}.centers_bkgd_ctlr_{level_tag}.tsv'), 
                          sep = "\t" )
        
        # Many of meta-clonotypes contain redundant information. 
        # We can winnow down to less-redundant list. We do this 
        # by ranking clonotypes from most to least specific. 
            # <min_nsubject> is minimum publicity of the meta-clonotype,  
            # <min_nr> is minimum non-redundancy
        # Add neighbors, K_neighbors, and nsubject columns
        from tcrdist.public import _neighbors_variable_radius, _neighbors_sparse_variable_radius
        centers_df['neighbors'] = _neighbors_variable_radius(pwmat=tr.pw_beta, 
                                                            radius_list = tr.clone_df['radius'])
        centers_df['K_neighbors'] = centers_df['neighbors'].apply(lambda x : len(x))
        # We determine how many <nsubjects> are in the set of neighbors 
        centers_df['nsubject']  = centers_df['neighbors'].\
                apply(lambda x: tr.clone_df['subject'].iloc[x].nunique())
        centers_df.to_csv(os.path.join(out_path, 
                                       f'{antigen_enriched_name}_{user_chain}.centers_bkgd_ctlr_{level_tag}.tsv'), 
                          sep = "\t" )

        from tcrdist.centers import rank_centers
        ranked_centers_df = rank_centers(
            centers_df = centers_df, 
            rank_column = 'chi2joint', 
            min_nsubject = min_subjects, 
            min_nr = min_redundancy)
        ############################################################################
        # Step 6.3:  (theta = 1E5) NR meta-clonotypes .tsv file                  ###
        ############################################################################
        # Output, ready to search bulk data.
        ranked_centers_df.to_csv( os.path.join(out_path, f'{antigen_enriched_name}_{user_chain}.ranked_centers_bkgd_ctlr_{level_tag}.tsv'), sep = "\t" )
        ############################################################################
        # Step 6.4: (theta = 1E5) Output Meta-Clonotypes HTML Summary            ###
        ############################################################################
        # Here we can make a svg logo for each NR meta-clonotype
        if ranked_centers_df.shape[0] > 0:
            from progress.bar import IncrementalBar
            from tcrdist.public import make_motif_logo
            svgs = list()
            svgs_raw = list()
            bar = IncrementalBar('Processing', max = ranked_centers_df.shape[0])
            for i,r in ranked_centers_df.iterrows():
                bar.next()
                centroid = r["cdr3_b_aa"]
                v_gene   = r["v_b_gene"]
                svg, svg_raw = make_motif_logo( tcrsampler = ts, 
                                                pwmat = tr.pw_beta,
                                                clone_df = tr.clone_df,
                                                centroid = centroid ,
                                                v_gene = v_gene ,
                                                radius = r['radius'],
                                                pwmat_str = 'pw_beta',
                                                cdr3_name = 'cdr3_b_aa',
                                                v_name = "v_b_gene",
                                                gene_names = ['v_b_gene','j_b_gene'])
                svgs.append(svg)
                svgs_raw.append(svg_raw)
            bar.next();bar.finish()
            ranked_centers_df['svg']      = svgs
            ranked_centers_df['svg_raw'] = svgs_raw

        def shrink(s):
            return s.replace('height="100%"', 'height="20%"').replace('width="100%"', 'width="20%"')
        labels =['cdr3_b_aa','v_b_gene', 'j_b_gene', 'pgen',
                'radius', 'regex','nsubject','K_neighbors', 
                'bkgd_hits_weighted','chi2dist','chi2re','chi2joint']
        
        output_html_name = os.path.join(out_path, 
                                        f'{antigen_enriched_name}_{user_chain}.ranked_centers_bkgd_ctlr_{level_tag}.html')

        with open(output_html_name, 'w') as output_handle:
            for i,r in ranked_centers_df.iterrows():
                #import pdb; pdb.set_trace()
                svg, svg_raw = r['svg'],r['svg_raw']
                output_handle.write("<br></br>")
                output_handle.write(shrink(svg))
                output_handle.write(shrink(svg_raw))
                output_handle.write("<br></br>")
                output_handle.write(pd.DataFrame(r[labels]).transpose().to_html())
                output_handle.write("<br></br>")
        # To ascertain which meta-clonotypes are likely to be most specific, 
        # take advantage of an existing function <bkgd_cntrl_nn2>.       

        ############################################################################
        # Step 6.5: Find optimal radii  (theta = 1E6)                            ###
        ############################################################################
        # level_tag = '1E6'
        # from tcrdist.neighbors import bkgd_cntl_nn2
        # centers_df  = bkgd_cntl_nn2(
        #     tr               = tr,
        #     tr_background    = tr_bkgd,
        #     weights          = tr_bkgd.clone_df.weights,
        #     ctrl_bkgd        = 10**-6, 
        #     col              = 'cdr3_b_aa',
        #     add_cols         = ['v_b_gene', 'j_b_gene'],
        #     ncpus            = 4,
        #     include_seq_info = True,
        #     thresholds       = [x for x in range(0,50,2)],
        #     generate_regex   = True,
        #     test_regex       = True,
        #     forced_max_radius = 36)
        # ############################################################################
        # # Step 6.6: (theta = 1E6) ALL meta-clonotypes .tsv file                   ##
        # ############################################################################
        # # save center to project_path for future use
        # centers_df.to_csv(os.path.join(out_path, 
        #                                f'{antigen_enriched_name}_{user_chain}.centers_bkgd_ctlr_{level_tag}.tsv'), 
        #                   sep = "\t" )
        
        # # Many of meta-clonotypes contain redundant information. 
        # # We can winnow down to less-redundant list. We do this 
        # # by ranking clonotypes from most to least specific. 
        #     # <min_nsubject> is minimum publicity of the meta-clonotype,  
        #     # <min_nr> is minimum non-redundancy
        # # Add neighbors, K_neighbors, and nsubject columns
        # from tcrdist.public import _neighbors_variable_radius, _neighbors_sparse_variable_radius
        # centers_df['neighbors'] = _neighbors_variable_radius(pwmat=tr.pw_beta, 
        #                                                      radius_list = centers_df['radius'])
        # centers_df['K_neighbors'] = centers_df['neighbors'].apply(lambda x : len(x))
        # # We determine how many <nsubjects> are in the set of neighbors 
        # centers_df['nsubject']  = centers_df['neighbors'].\
        #     apply(lambda x: tr.clone_df['subject'].iloc[x].nunique())
        # centers_df.to_csv(os.path.join(out_path, 
        #                                f'{antigen_enriched_name}_{user_chain}.centers_bkgd_ctlr_{level_tag}.tsv'), 
        #                   sep = "\t" )

        # from tcrdist.centers import rank_centers
        # ranked_centers_df = rank_centers(
        #     centers_df = centers_df, 
        #     rank_column = 'chi2joint', 
        #     min_nsubject = min_subjects, 
        #     min_nr = min_redundancy)
        # ############################################################################
        # # Step 6.7:  (theta = 1E6) NR meta-clonotypes .tsv file                  ###
        # ############################################################################
        # # Output, ready to search bulk data.
        # ranked_centers_df.to_csv(os.path.join(out_path, 
        #                                       f'{antigen_enriched_name}_{user_chain}.ranked_centers_bkgd_ctlr_{level_tag}.tsv'), 
        #                          sep = "\t")

        # ############################################################################
        # # Step 6.8: (theta = 1E6) Output Meta-Clonotypes HTML Summary            ###
        # ############################################################################
        # # Here we can make a svg logo for each meta-clonotype
        # from progress.bar import IncrementalBar
        # from tcrdist.public import make_motif_logo
        # if ranked_centers_df.shape[0] > 0:
        #     svgs = list()
        #     svgs_raw = list()
        #     bar = IncrementalBar('Processing', max = ranked_centers_df.shape[0])
        #     for i,r in ranked_centers_df.iterrows():
        #         bar.next()
        #         centroid = r["cdr3_b_aa"]
        #         v_gene   = r["v_b_gene"]
        #         svg, svg_raw = make_motif_logo( tcrsampler = ts, 
        #                                         pwmat = tr.pw_beta,
        #                                         clone_df = tr.clone_df,
        #                                         centroid = centroid ,
        #                                         v_gene = v_gene ,
        #                                         radius = r['radius'],
        #                                         pwmat_str = 'pw_beta',
        #                                         cdr3_name = 'cdr3_b_aa',
        #                                         v_name = 'v_b_gene',
        #                                         gene_names = ['v_b_gene','j_b_gene'])
        #         svgs.append(svg)
        #         svgs_raw.append(svg_raw)
        #     bar.next();bar.finish()
        #     ranked_centers_df['svg']      = svgs
        #     ranked_centers_df['svg_raw'] = svgs_raw

        #     def shrink(s):
        #         return s.replace('height="100%"', 'height="20%"').replace('width="100%"', 'width="20%"')
        #     labels =['cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'pgen', 'radius', 'regex','nsubject','K_neighbors', 'bkgd_hits_weighted','chi2dist','chi2re','chi2joint']
            
        #     output_html_name = os.path.join(out_path, f'{antigen_enriched_name}_{user_chain}.ranked_centers_bkgd_ctlr_{level_tag}.html')

        #     with open(output_html_name, 'w') as output_handle:
        #         for i,r in ranked_centers_df.iterrows():
        #             #import pdb; pdb.set_trace()
        #             svg, svg_raw = r['svg'],r['svg_raw']
        #             output_handle.write("<br></br>")
        #             output_handle.write(shrink(svg))
        #             output_handle.write(shrink(svg_raw))
        #             output_handle.write("<br></br>")
        #             output_handle.write(pd.DataFrame(r[labels]).transpose().to_html())
        #             output_handle.write("<br></br>")







    #### Running code for alpha chain below






    elif user_chain == "alpha":
        print(f"USING tcrsampler TO CONSTRUCT A CUSTOM V-J MATCHED BACKGROUND")
        from tcrsampler.sampler import TCRsampler
        ## I had to change the default_background of for the alpha chains. 
        ## The backgrounds are supplied by TCRsampler setup_db.py script. 
        ## They are stored wherever the tcrsampler package is installed. 
        ## For John, that is at /stor/home/jfm2773/anaconda3/envs/tcrdist3_code/lib/python3.10/site-packages/tcrsampler/db
        ts = TCRsampler(default_background = 'ruggiero_human_alpha_t.tsv.sampler.tsv')
        ## This reference contains alpha chains from healthy patients from this paper: doi: 10.1038/ncomms9081 
        ## The healthy donors are age 26 - 65. 4 females, 2 males, 2 HDs are CMV+, some HLA diversity, but only 6 patients so still very limited.
        ## Important to note that there aren't any kids to be used for comparison in this dataset. 
        # Stratify sample so that each subject contributes similarly to estimate of 
        # gene usage frequency
        from tcrdist.background import get_stratified_gene_usage_frequency
        ts = get_stratified_gene_usage_frequency(ts = ts, replace = True) 
        # Synthesize an inverse probability weighted V,J gene background that matches 
        # usage in your enriched repertoire 
        df_vj_background = tr.synthesize_vj_matched_background(ts = ts, chain = user_chain)
        
        df_chosen_bkgd = sample_df_bkgd(ts= ts, 
                                        nsubject= 6,
                                        chain= user_chain,
                                        size= 9000)
            ## This is my general version of sample_britanova().
            ## Pass it ts from ts= get_stratified_gene_usage_frequency() w/your background of choice.
        
        df_chosen_bkgd = get_gene_frequencies(ts= ts, 
                                              df= df_chosen_bkgd,
                                              cols= ["v_a_gene", "j_a_gene"])
        df_chosen_bkgd['weights'] = 1
        df_chosen_bkgd['source']  = "stratified_random"
        #print(df_chosen_bkgd.head())
    
        # Combine the two parts of the background into a single DataFrame
        df_bkgd= pd.concat([df_vj_background.copy(),
                           df_chosen_bkgd.copy()]).\
                           reset_index(drop= True)

        # Save the background for future use
        background_outfile = os.path.join(out_path, 
                                          f"{antigen_enriched_name}_{user_chain}.olga100K_bkgd.csv")
        print(f'WRITING {background_outfile}')
        df_bkgd.to_csv(background_outfile, index = False)
        # Load the background to a TCRrep without computing pairwise distances 
        # (i.e., compute_distances = False)
        tr_bkgd = TCRrep(
            cell_df = df_bkgd,
            organism = "human", 
            chains = [user_chain], 
            compute_distances = False)
    
        ############################################################################
        # Step 4: Calculate Distances                                          #####
        ############################################################################
        print(f"COMPUTING RECTANGULARE DISTANCE")

        tr.compute_sparse_rect_distances(
            df = tr.clone_df,
            df2 = tr_bkgd.clone_df,
            radius=50,
            chunk_size = 100)
            ## This is the function that takes a while. 
        scipy_path= os.path.join(out_path, 
                                f"{antigen_enriched_name}_{user_chain}.rw_alpha.npz") 
        scipy.sparse.save_npz(scipy_path, tr.rw_alpha)
            # Tip: For larger dataset you can use a sparse implementation: 
            # 30.8 s ± 0 ns per loop ; tr.cpus = 6
            # %timeit -r tr.compute_sparse_rect_distances(df = tr.clone_df, df2 = tr_bkdg.clone_df,radius=50, chunk_size=85)
        ############################################################################
        # Step 5: Examine Density ECDFS                                        #####
        ############################################################################
            # Investigate the density of neighbors to each TCR, based on expanding 
            # distance radius.
        from tcrdist.ecdf import distance_ecdf, _plot_manuscript_ecdfs
        import matplotlib.pyplot as plt
            # Compute empirical cumulative density function (ecdf)
            # Compare Antigen Enriched TCRs (against itself).
        thresholds, antigen_enriched_ecdf = distance_ecdf(
            tr.pw_alpha,    
            thresholds=range(0,50,2))
            # Compute empirical cumulative density function (ecdf)
            # Compare Antigen Enriched TCRs (against) 200K probability 
            # inverse weighted background
        thresholds, background_ecdf = distance_ecdf(
            tr.rw_alpha,
            thresholds=range(0,50,2),
            weights= tr_bkgd.clone_df['weights'], 
            absolute_weight = True)
            # plot_ecdf similar to tcrdist3 manuscript #
        antigen_enriched_ecdf[antigen_enriched_ecdf == antigen_enriched_ecdf.min()] = 1E-10
        f1 = _plot_manuscript_ecdfs(
            thresholds, 
            antigen_enriched_ecdf, 
            ylab= 'Proportion of Antigen Enriched TCRs', 
            cdr3_len=tr.clone_df.cdr3_a_aa.str.len(), 
            min_freq=1E-10)
        f1.savefig(os.path.join(out_path, f'{antigen_enriched_name}_{user_chain}.ecdf_AER_plot.png'))
        f2 = _plot_manuscript_ecdfs(
            thresholds,
            background_ecdf,
            ylab= 'Proportion of Reference TCRs',
            cdr3_len=tr.clone_df.cdr3_a_aa.str.len(),
            min_freq=1E-10)
        f2.savefig(os.path.join(out_path, f'{antigen_enriched_name}_{user_chain}.ecdf_BUR_plot.png'))

        ############################################################################
        # Step 6: Find optimal radii  (theta = 1E5                             #####
        ############################################################################
        # To ascertain which meta-clonotypes are likely to be most specific, 
        # take advantage of an existing function <bkgd_cntrl_nn2>.                                                                                                                                  

        level_tag = '1E5'
        from tcrdist.neighbors import bkgd_cntl_nn2
        centers_df  = bkgd_cntl_nn2(
            tr               = tr,
            tr_background    = tr_bkgd,
            weights          = tr_bkgd.clone_df.weights,
            ctrl_bkgd        = 10**-5, 
            col              = 'cdr3_a_aa',
            add_cols         = ['v_a_gene', 'j_a_gene'],
            pw_mat_str       = 'pw_alpha', 
            rw_mat_str       = 'rw_alpha',
            ncpus            = 4,
            include_seq_info = True,
            thresholds       = [x for x in range(0,50,2)],
            generate_regex   = True,
            test_regex       = True,
            forced_max_radius = 36)

        ############################################################################
        # Step 6.2: (theta = 1E5) ALL meta-clonotypes .tsv file                   ##
        ############################################################################
        # save center to project_path for future use
        centers_df.to_csv(os.path.join(out_path, 
                                       f'{antigen_enriched_name}.centers_bkgd_ctlr_{level_tag}.tsv'), 
                          sep = "\t" )
        
        # Many of meta-clonotypes contain redundant information. 
        # We can winnow down to less-redundant list. We do this 
        # by ranking clonotypes from most to least specific. 
            # <min_nsubject> is minimum publicity of the meta-clonotype,  
            # <min_nr> is minimum non-redundancy
        # Add neighbors, K_neighbors, and nsubject columns
        from tcrdist.public import _neighbors_variable_radius, _neighbors_sparse_variable_radius
        centers_df['neighbors'] = _neighbors_variable_radius(pwmat=tr.pw_alpha, 
                                                            radius_list = tr.clone_df['radius'])
        centers_df['K_neighbors'] = centers_df['neighbors'].apply(lambda x : len(x))
        # We determine how many <nsubjects> are in the set of neighbors 
        centers_df['nsubject']  = centers_df['neighbors'].\
                apply(lambda x: tr.clone_df['subject'].iloc[x].nunique())
        centers_df.to_csv(os.path.join(out_path, 
                                       f'{antigen_enriched_name}_{user_chain}.centers_bkgd_ctlr_{level_tag}.tsv'), 
                          sep = "\t" )

        from tcrdist.centers import rank_centers
        ranked_centers_df = rank_centers(
            centers_df = centers_df, 
            rank_column = 'chi2joint', 
            min_nsubject = min_subjects, 
            min_nr = min_redundancy)
        ############################################################################
        # Step 6.3:  (theta = 1E5) NR meta-clonotypes .tsv file                  ###
        ############################################################################
        # Output, ready to search bulk data.
        ranked_centers_df.to_csv(os.path.join(out_path, 
                                              f'{antigen_enriched_name}_{user_chain}.ranked_centers_bkgd_ctlr_{level_tag}.tsv'), 
                                 sep = "\t" )
        ############################################################################
        # Step 6.4: (theta = 1E5) Output Meta-Clonotypes HTML Summary            ###
        ############################################################################
        # Here we can make a svg logo for each NR meta-clonotype
        if ranked_centers_df.shape[0] > 0:
            from progress.bar import IncrementalBar
            from tcrdist.public import make_motif_logo
            svgs = list()
            svgs_raw = list()
            bar = IncrementalBar('Processing', max = ranked_centers_df.shape[0])
            for i,r in ranked_centers_df.iterrows():
                bar.next()
                centroid = r["cdr3_a_aa"]
                v_gene   = r["v_a_gene"]
                svg, svg_raw = make_motif_logo( tcrsampler = ts, 
                                                pwmat = tr.rw_alpha,
                                                clone_df = tr.clone_df,
                                                centroid = centroid ,
                                                v_gene = v_gene ,
                                                radius = r['radius'],
                                                pwmat_str = 'rw_alpha',
                                                cdr3_name = 'cdr3_a_aa',
                                                v_name = "v_a_gene",
                                                gene_names = ['v_a_gene','j_a_gene'])
                svgs.append(svg)
                svgs_raw.append(svg_raw)
            bar.next();bar.finish()
            ranked_centers_df['svg']      = svgs
            ranked_centers_df['svg_raw'] = svgs_raw

        def shrink(s):
            return s.replace('height="100%"', 'height="20%"').replace('width="100%"', 'width="20%"')
        labels =['cdr3_a_aa','v_a_gene', 'j_a_gene', 'pgen',
                'radius', 'regex','nsubject','K_neighbors', 
                'bkgd_hits_weighted','chi2dist','chi2re','chi2joint']
        
        output_html_name = os.path.join(out_path, 
                                        f'{antigen_enriched_name}_{user_chain}.ranked_centers_bkgd_ctlr_{level_tag}.html')

        with open(output_html_name, 'w') as output_handle:
            for i,r in ranked_centers_df.iterrows():
                #import pdb; pdb.set_trace()
                svg, svg_raw = r['svg'],r['svg_raw']
                output_handle.write("<br></br>")
                output_handle.write(shrink(svg))
                output_handle.write(shrink(svg_raw))
                output_handle.write("<br></br>")
                output_handle.write(pd.DataFrame(r[labels]).transpose().to_html())
                output_handle.write("<br></br>")
        # To ascertain which meta-clonotypes are likely to be most specific, 
        # take advantage of an existing function <bkgd_cntrl_nn2>.       

        # ############################################################################
        # # Step 6.5: Find optimal radii  (theta = 1E6)                            ###
        # ############################################################################
        # level_tag = '1E6'
        # from tcrdist.neighbors import bkgd_cntl_nn2
        # centers_df  = bkgd_cntl_nn2(
        #     tr               = tr,
        #     tr_background    = tr_bkgd,
        #     weights          = tr_bkgd.clone_df.weights,
        #     ctrl_bkgd        = 10**-6, 
        #     col              = 'cdr3_a_aa',
        #     add_cols         = ['v_a_gene', 'j_a_gene'],
        #     pw_mat_str       = 'pw_alpha', 
        #     rw_mat_str       = 'rw_alpha',
        #     ncpus            = 4,
        #     include_seq_info = True,
        #     thresholds       = [x for x in range(0,50,2)],
        #     generate_regex   = True,
        #     test_regex       = True,
        #     forced_max_radius = 36)
        # ############################################################################
        # # Step 6.6: (theta = 1E6) ALL meta-clonotypes .tsv file                   ##
        # ############################################################################
        # # save center to project_path for future use
        # centers_df.to_csv(os.path.join(out_path, 
        #                                f'{antigen_enriched_name}_{user_chain}.centers_bkgd_ctlr_{level_tag}.tsv'), 
        #                   sep = "\t" )
        
        # # Many of meta-clonotypes contain redundant information. 
        # # We can winnow down to less-redundant list. We do this 
        # # by ranking clonotypes from most to least specific. 
        #     # <min_nsubject> is minimum publicity of the meta-clonotype,  
        #     # <min_nr> is minimum non-redundancy
        # # Add neighbors, K_neighbors, and nsubject columns
        # from tcrdist.public import _neighbors_variable_radius, _neighbors_sparse_variable_radius
        # centers_df['neighbors'] = _neighbors_variable_radius(pwmat=tr.pw_alpha, 
        #                                                      radius_list = centers_df['radius'])
        # centers_df['K_neighbors'] = centers_df['neighbors'].apply(lambda x : len(x))
        # # We determine how many <nsubjects> are in the set of neighbors 
        # centers_df['nsubject']  = centers_df['neighbors'].\
        #     apply(lambda x: tr.clone_df['subject'].iloc[x].nunique())
        # centers_df.to_csv(os.path.join(out_path, 
        #                                f'{antigen_enriched_name}_{user_chain}.centers_bkgd_ctlr_{level_tag}.tsv'), 
        #                   sep = "\t" )

        # from tcrdist.centers import rank_centers
        # ranked_centers_df = rank_centers(
        #     centers_df = centers_df, 
        #     rank_column = 'chi2joint', 
        #     min_nsubject = min_subjects, 
        #     min_nr = min_redundancy)
        # ############################################################################
        # # Step 6.7:  (theta = 1E6) NR meta-clonotypes .tsv file                  ###
        # ############################################################################
        # # Output, ready to search bulk data.
        # ranked_centers_df.to_csv(os.path.join(out_path, 
        #                                       f'{antigen_enriched_name}_{user_chain}.ranked_centers_bkgd_ctlr_{level_tag}.tsv'), 
        #                          sep = "\t")

        # ############################################################################
        # # Step 6.8: (theta = 1E6) Output Meta-Clonotypes HTML Summary            ###
        # ############################################################################
        # # Here we can make a svg logo for each meta-clonotype
        # from progress.bar import IncrementalBar
        # from tcrdist.public import make_motif_logo
        # if ranked_centers_df.shape[0] > 0:
        #     cdr3_name = 'cdr3_a_aa'
        #     svgs = list()
        #     svgs_raw = list()
        #     bar = IncrementalBar('Processing', max = ranked_centers_df.shape[0])
        #     for i,r in ranked_centers_df.iterrows():
        #         bar.next()
        #         centroid = r[cdr3_name]
        #         v_gene   = r['v_a_gene']
        #         svg, svg_raw = make_motif_logo( tcrsampler = ts, 
        #                                         pwmat = tr.pw_alpha,
        #                                         clone_df = tr.clone_df,
        #                                         centroid = centroid ,
        #                                         v_gene = v_gene ,
        #                                         radius = r['radius'],
        #                                         pwmat_str = 'pw_alpha',
        #                                         cdr3_name = 'cdr3_a_aa',
        #                                         v_name = 'v_a_gene',
        #                                         gene_names = ['v_a_gene','j_a_gene'])
        #         svgs.append(svg)
        #         svgs_raw.append(svg_raw)
        #     bar.next();bar.finish()
        #     ranked_centers_df['svg']      = svgs
        #     ranked_centers_df['svg_raw'] = svgs_raw

        #     def shrink(s):
        #         return s.replace('height="100%"', 'height="20%"').replace('width="100%"', 'width="20%"')
        #     labels =['cdr3_a_aa', 'v_a_gene', 'j_a_gene', 'pgen', 'radius', 'regex','nsubject','K_neighbors', 'bkgd_hits_weighted','chi2dist','chi2re','chi2joint']
            
        #     output_html_name = os.path.join(out_path, 
        #                                     f'{antigen_enriched_name}_{user_chain}.ranked_centers_bkgd_ctlr_{level_tag}.html')

        #     with open(output_html_name, 'w') as output_handle:
        #         for i,r in ranked_centers_df.iterrows():
        #             #import pdb; pdb.set_trace()
        #             svg, svg_raw = r['svg'],r['svg_raw']
        #             output_handle.write("<br></br>")
        #             output_handle.write(shrink(svg))
        #             output_handle.write(shrink(svg_raw))
        #             output_handle.write("<br></br>")
        #             output_handle.write(pd.DataFrame(r[labels]).transpose().to_html())
        #             output_handle.write("<br></br>")
            ## I didn't fined 1E6 runs useful, so to decrease computation, I've commented these out. 
    else:
        raise ValueError("Unknown user_chain parameter. Only 'alpha' and 'beta' are accepted.")
    


if __name__ == "__main__":
  
    find_metaclonotypes(user_chain= "alpha",
                        input_file= '/stor/work/Ehrlich_COVID19/Users/John/projects/stim_TCRs/TCRdist/tcrdist3_analysis/input_data/stim_and_ref_TCRs_072222.csv',
                        out_path= '/stor/work/Ehrlich_COVID19/Users/John/projects/stim_TCRs/TCRdist/tcrdist3_analysis/output/stimRef_sub2_red1',
                        ncpus= 10,
                        min_subjects= 2,
                        min_redundancy= 1)
    
    find_metaclonotypes(user_chain= "beta",
                        input_file= '/stor/work/Ehrlich_COVID19/Users/John/projects/stim_TCRs/TCRdist/tcrdist3_analysis/input_data/stim_and_ref_TCRs_072222.csv',
                        out_path= '/stor/work/Ehrlich_COVID19/Users/John/projects/stim_TCRs/TCRdist/tcrdist3_analysis/output/stimRef_sub2_red1',
                        ncpus= 10,
                        min_subjects= 2,
                        min_redundancy= 1)
