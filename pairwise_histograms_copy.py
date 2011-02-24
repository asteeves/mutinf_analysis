## Program to plot 2D histograms for torsions with highest mutual information given a pair of residues

import sys
sys.path.append('/home/mcclendon/trajtools/')
from dihedral_mutent import *
import matplotlib
import matplotlib.pyplot as plt

def find_in_resfile_and_append(reslist, run_params, res_name_i, res_num_i):
        rp = run_params
        resfile=open(rp.resfile_fn,'r')
        reslines=resfile.readlines()
        resfile.close()
        all_angle_info=None
        for resline in reslines:
            if len(resline.strip()) == 0: continue
            xvg_resnum, res_name, res_num = resline.split()
            if res_num == res_num_i and res_name == res_name_i: 
                reslist.append(ResidueChis(res_name,res_num, xvg_resnum, rp.xvg_basedir, rp.num_sims, rp.num_structs, rp.xvgorpdb, rp.binwidth, rp.sigalpha, rp.permutations, rp.phipsi, rp.backbone_only, rp.adaptive_partitioning, rp.which_runs, rp.pair_runs, bootstrap_choose = rp.bootstrap_choose, calc_variance=rp.calc_variance, all_angle_info=all_angle_info, xvg_chidir=rp.xvg_chidir, skip=rp.skip,skip_over_steps=rp.skip_over_steps, calc_mutinf_between_sims=rp.calc_mutinf_between_sims,max_num_chis=rp.max_num_chis))
        return reslist

class Mutinf_Matrix_Chis_Bootstraps:
    mut_info_res_matrix = None
    rownames = []
    colnames = []
    old_bootstrap_sets = 0
    prefix = ""

    def __init__(self, prefix):
        self.old_bootstrap_sets = 0 #figure out # of bootstrap sets originally from output data
        self.prefix = prefix
        while (os.path.exists(self.prefix+"_bootstrap_"+str(self.old_bootstrap_sets)+"_mutinf_0diag.txt")):
           self.old_bootstrap_sets += 1
        if self.old_bootstrap_sets == 0:
           print "cannot find file"+self.prefix+"run1"+str(0)+"_mutinf_0diag.txt"
        print "old bootstrap sets: "+str(self.old_bootstrap_sets)
        ## READ MATRIX
        # self.prefix = "%s-nsims%d-structs%d-bin%d" % (os.path.basename(self.resfile_fn), self.num_sims, self.num_structs, int(self.binwidth))
        (test_matrix, self.rownames, self.colnames) = read_matrix_chis(self.prefix+"_bootstrap_0_mutinf_0diag.txt")
        self.rownames = self.colnames
        ##print test_matrix
        self.mut_info_res_matrix = zeros((self.old_bootstrap_sets, test_matrix.shape[0],test_matrix.shape[1],test_matrix.shape[2],test_matrix.shape[3]),float64)
        for bootstrap in range(self.old_bootstrap_sets):
            (self.mut_info_res_matrix[bootstrap,:,:,:,:], self.rownames, self.colnames) = read_matrix_chis(self.prefix+"_bootstrap_"+str(bootstrap)+"_mutinf_0diag.txt")

    def pairwise_histogram(self,res_i_index, res_j_index, xvg_basedir, xvg_chidir="/", skip=1):
       ### Setup options so we can use class ResidueChis from dihedral_mutent.py
       prefixline = re.compile(r'(.*).nsims([0-9]+).structs([0-9]+).bin([0-9]+)$')
       matches_prefixline = prefixline.match(self.prefix)
       print matches_prefixline
       resfile_fn = matches_prefixline.group(1)
       num_sims = int(matches_prefixline.group(2))
       num_structs = int(matches_prefixline.group(3))
       binwidth= int(matches_prefixline.group(4))
       sigalpha = 0.01
       permutations = 0
       adaptive = "no"
       backbone = "phipsichi"
       bootstrap_set_size = num_sims
       correct_formutinf_between_sims = "no"
       plot_2d_histograms = "no" #since we're doing it outside the mutinf code
       #compiler = 'gcc'
       adaptive_partitioning = 0
       phipsi = 0
       zoom_to_step=0
       permutation = 0
       bootstrap_sets = 1
       backbone_only = 0
       max_num_chis= 99
       sigalpha=0.01
       correct_formutinf_between_sims = "no"
       load_matrices_numstructs = 0

       if backbone == "phipsichi":
           phipsi = 2;
           backbone_only =0
       if backbone == "phipsi":
           phipsi = 2;
           backbone_only = 1
       #phipsi = options.backbone.find("phipsi")!=-1
       #backbone_only = options.backbone.find("chi")==-1
       bins=arange(-180,180,binwidth) #Compute bin edges
       nbins = len(bins)
       xvgorpdb = "xvg"
       traj_fns = None

       #### SETUP BOOTSTRAP RESAMPLING 
       if bootstrap_set_size == None:
           bootstrap_set_size = num_sims

       which_runs = []
       pair_runs_list = []
       for myruns in xuniqueCombinations(range(num_sims), bootstrap_set_size):
           which_runs.append(myruns)
       print which_runs

       bootstrap_pair_runs_list = []
       for bootstrap in range((array(which_runs)).shape[0]):
         pair_runs_list = []
         for myruns2 in xcombinations(which_runs[bootstrap], 2):
             pair_runs_list.append(myruns2)
         bootstrap_pair_runs_list.append(pair_runs_list)

       pair_runs_array = array(bootstrap_pair_runs_list,int16)
       print pair_runs_array 

       
       #### DONE WITH SETUP BOOTSTRAP RESAMPLING 


       ## SETUP RUN PARAMETERS
       run_params = RunParameters(resfile_fn=resfile_fn, adaptive_partitioning=adaptive_partitioning, phipsi=phipsi, backbone_only=backbone_only, nbins = nbins,
           bootstrap_set_size=bootstrap_set_size, sigalpha=sigalpha, permutations=permutations, num_sims=num_sims, num_structs=num_structs,
           binwidth=binwidth, bins=bins, which_runs=which_runs, xvgorpdb=xvgorpdb, traj_fns=traj_fns, xvg_basedir=xvg_basedir, calc_variance=False, xvg_chidir=xvg_chidir,bootstrap_choose=bootstrap_set_size,pair_runs=pair_runs_array,skip=skip,skip_over_steps=zoom_to_step,calc_mutinf_between_sims=correct_formutinf_between_sims,load_matrices_numstructs=load_matrices_numstructs,plot_2d_histograms=plot_2d_histograms,max_num_chis=max_num_chis)
       ## DONE SETTING UP RUN PARAMETERS

       
           
       #mutinf_boots = self.mut_info_res_matrix[:,i,j,k,m].copy()

       #### Average over non-zero mutinf values from bootstraps
       mutinf_pairs_these_residues = zeros((self.old_bootstrap_sets,6,6),float64)
       mutinf_pairs_these_residues_avg = zeros((6,6),float64)
       #print "self.mut_info_res_matrix shape:" + str(self.mut_info_res_matrix_avg.shape)
       print "residue i index: "+str(res_i_index)+" residue j index: "+str(res_j_index)
       mutinf_pairs_these_residues[:,:,:] = self.mut_info_res_matrix[:,int(res_i_index),int(res_j_index),:,:]
       print mutinf_pairs_these_residues

       mutinf_pairs_these_residues_avg_list = []
       for i in range(6):
           for j in range(6):
               if(sum(mutinf_pairs_these_residues[:,i,j]) > 0):
                   mutinf_boots = mutinf_pairs_these_residues[:,i,j].copy()
                   mutinf_boots[mutinf_boots < 0] = 0 #use negative and zero values for significance testing(if desired) but not in the average
                   #zero values of mutinf_boots will include those zeroed out in the permutation test
                   mutinf_pairs_these_residues_avg[i,j] = average(mutinf_boots[mutinf_boots > 0 ],axis=0)
                   mutinf_pairs_these_residues_avg_list.append(mutinf_pairs_these_residues_avg[i,j])
               else:
                   mutinf_pairs_these_residues_avg_list.append(0.0)
       maxchi = 6
       print self.rownames
       print self.colnames
       res_i_line = self.rownames[int(res_i_index)]  #because there are six entries for each residue
       res_j_line = self.colnames[int(res_j_index)]  #because there are six entries for each residue

       mutinf_pairs_these_residues_avg_list.sort()
       mutinf_pairs_these_residues_avg_list.reverse()
       ## find top three mutinf pairs in matrix
       sorted_mutinfs = mutinf_pairs_these_residues_avg_list
       print "sorted mutinfs: " + str(sorted_mutinfs)
       top_mutinfs = sorted_mutinfs[0:3]
       print "top mutinfs: " + str(top_mutinfs)
       top_mutinfs_chi_i = zeros((3),float64)
       top_mutinfs_chi_j = zeros((3),float64)
       print "mutinf pairs submatrix:"
       print mutinf_pairs_these_residues_avg
       counter = 0
       mutinf_pairs_these_residues_avg_flat = (mutinf_pairs_these_residues_avg.copy()).flatten()       
       for counter in range(3):
           print nonzero(mutinf_pairs_these_residues_avg_flat == top_mutinfs[counter])
           ij_top_mutinf = nonzero(mutinf_pairs_these_residues_avg_flat == top_mutinfs[counter])
           print "chi index:"+str(ij_top_mutinf)
           top_mutinfs_chi_i[counter] = int(ij_top_mutinf[0]) / 6
           top_mutinfs_chi_j[counter] = int(ij_top_mutinf[0]) % 6
       
       #for i in range(6):
       #    for j in range(6):
       #        if counter <=2: 
       #            if mutinf_pairs_these_residues_avg[i,j] < top_mutinfs[counter]+0.00001 and  mutinf_pairs_these_residues_avg[i,j] >  top_mutinfs[counter]-0.00001:
       #                print "mutinf value:"+str(top_mutinfs[counter])
       #                print "chi i:       "+str(i)
       #                print "chi j:       "+str(j)
       #                top_mutinfs_chi_i[counter] = i
       #                top_mutinfs_chi_j[counter] = j
       #                counter += 1
       
       #### Grab data and make marginal distributions for selected residues using class ResidueChis from dihedral_mutent.py

       resline_regexp = re.compile(r'([A-Z]+)([0-9]+[A-Z]*)') #'_([0-9]+)')
       
       print "residue: "+res_i_line
       matches_i = resline_regexp.match(res_i_line)
       res_i_type = matches_i.group(1)
       res_i_num = matches_i.group(2) #includes number followed by chain name
       
       matches_j = resline_regexp.match(res_j_line)
       res_j_type = matches_j.group(1)
       res_j_num = matches_j.group(2)
       #res_j_chain = matches_j.group(2)
       
       #pdb_lines = commands.getoutput("egrep '^ATOM ' %s" % options.structure).split("\n")
       #pdb_obj = PDBlite.PDB(pdb_lines,fn=options.structure)
       #print pdb_obj
       
       reslist = []
       find_in_resfile_and_append(reslist, run_params, res_i_type, res_i_num)
       find_in_resfile_and_append(reslist, run_params, res_j_type, res_j_num)
       
       #### Done grabbing data
       
       Pij_top_mutinfs = []
       PiPj_top_mutinfs = []
       for counter_top_mutinfs in range(3):
           mychi1 = top_mutinfs_chi_i[counter_top_mutinfs]
           mychi2 = top_mutinfs_chi_j[counter_top_mutinfs]
           count_matrix = zeros((bootstrap_sets, permutations + 1 , nbins*nbins), int32)
           bins1 = reslist[0].bins[mychi1,:,:,:]
           bins2 = reslist[1].bins[mychi2,:,:,:]
           max_num_angles = int(max(reslist[0].numangles))
           bootstrap_choose = num_sims
           numangles_bootstrap = zeros((bootstrap_sets),int64)
           numangles_bootstrap[:] = int(sum(reslist[0].numangles))
           
           code = """
           // weave6 in pairwise histograms
           // bins dimensions: (permutations + 1) * bootstrap_sets * bootstrap_choose * max_num_angles
            //#include <math.h>
            for(int mybootstrap=0; mybootstrap < bootstrap_sets; mybootstrap++) {
             int mynumangles = 0;
             for (int permut=0; permut < permutations + 1; permut++) {
                 mynumangles = *(numangles_bootstrap + mybootstrap);
                 for (int anglenum=0; anglenum< mynumangles; anglenum++) {
                 if(mybootstrap == bootstrap_sets - 1) {
                   //printf("bin12 %i \\n",(*(bins1  +  mybootstrap*bootstrap_choose*max_num_angles  +  anglenum))*nbins +   (*(bins2 + permut*bootstrap_sets*bootstrap_choose*max_num_angles + mybootstrap*bootstrap_choose*max_num_angles  +  anglenum)));
                   }
                    *(count_matrix  +  mybootstrap*(permutations + 1)*nbins*nbins  +  permut*nbins*nbins  +  (*(bins1  +  mybootstrap*bootstrap_choose*max_num_angles  +  anglenum))*nbins +   (*(bins2 + permut*bootstrap_sets*bootstrap_choose*max_num_angles + mybootstrap*bootstrap_choose*max_num_angles  +  anglenum))) += 1;

                 }
               }
              }
             """
           if(VERBOSE >= 2): print "about to populate count_matrix"
           weave.inline(code, ['num_sims', 'numangles_bootstrap', 'nbins', 'bins1', 'bins2', 'count_matrix','bootstrap_sets','permutations','max_num_angles','bootstrap_choose'],
                    #type_converters = converters.blitz,
                    compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])
           
           #consider using first row and column for marginal distributions?
           Pij, PiPj = zeros((nbins+1, nbins+1), float64) - 1, zeros((nbins+1, nbins+1), float64) - 1
           Pij[1:,1:]  = (count_matrix[0,permutation,:]).reshape((nbins,nbins))
           marginals_i = sum(Pij[1:,1:],axis=1)
           marginals_j = sum(Pij[1:,1:],axis=0)
           PiPj[1:,1:] = outer(marginals_i + 0.0 , marginals_j + 0.0)
           #PiPj[1:,1:] = (ninj_flat[0,permutation,:]).reshape((nbins,nbins)) / (numangles_bootstrap[0] * 1.0)
           Pij_top_mutinfs.append(Pij)
           PiPj_top_mutinfs.append(PiPj)
       return Pij_top_mutinfs, PiPj_top_mutinfs, top_mutinfs_chi_i, top_mutinfs_chi_j, top_mutinfs, bins, self.rownames, self.colnames
   
   
    def plot_pairwise_histogram(self,matrix,bins,rownames,colnames):
       fig = plt.figure(figsize=(8,10),facecolor='w')
       ax1 = fig.add_axes([0.1,0.85,0.8,0.08],frameon=False)
       #dend = sch.dendrogram(links,orientation='top')
       ax1.set_xticks([])
       ax1.set_yticks([])
       axmatrix = fig.add_axes([0.1,0.1,0.8,0.7])
       #idx = dend['leaves']
       #reordered = self.mymatrix[idx,:]
       #reordered = reordered[:,idx]
       im = axmatrix.matshow(matrix, aspect='auto', origin='upper', cmap=plt.cm.binary)
       #im.set_clim(0.0,max(matrix[:,:]))
       axmatrix.set_xticks([])
       axmatrix.set_xticklabels(arange(0,360,360/len(bins)),rotation='vertical')
       axmatrix.set_yticks([])
       axmatrix.set_yticklabels(arange(0,360,360/len(bins)),rotation='vertical')
       #for tick in axmatrix.xaxis.get_major_ticks():
           #tick.tick1On = False
           #tick.tick2On = False
       for label in axmatrix.xaxis.get_ticklabels():
           label.set_rotation(90)
           label.set_size(6)
       for label in axmatrix.yaxis.get_ticklabels():
           label.set_rotation(90)
           label.set_size(6)
       fig.suptitle('Pairwise Histogram',size=12.0)
       return fig

if __name__ == "__main__":
    usage="%prog resfile residue_i residue_j  # where resfile is in the format <1-based-index> <aa type> <res num>"
    parser=OptionParser(usage)
    parser.add_option("-x", "--xvg_basedir", default=None, type="string", help="basedir to look for xvg files")
    #parser.add_option("-s", "--sigalpha", default=0.01, type="float", help="p-value threshold for statistical filtering, lower is stricter")
    parser.add_option("-s", "--num_structs", default=10001, type=int, help="number of snapshots per sim")
    parser.add_option("-w", "--binwidth", default=15.0, type="float", help="width of the bins in degrees")
    parser.add_option("-n", "--num_sims", default=None, type="int", help="number of simulations")
    #parser.add_option("-p", "--permutations", default=0, type="int", help="number of permutations for independent mutual information, for subtraction from total Mutual Information")
    parser.add_option("-d", "--xvg_chidir", default = "/dihedrals/g_chi/", type ="string", help="subdirectory under xvg_basedir/run# where chi angles are stored")
    #parser.add_option("-a", "--adaptive", default = "yes", type ="string", help="adaptive partitioning (yes|no)")
    #parser.add_option("-b", "--backbone", default = "phipsichi", type = "string", help="chi: just sc  phipsi: just bb  phipsichi: bb + sc")
    #parser.add_option("-o", "--bootstrap_set_size", default = None, type = "int", help="perform bootstrapping within this script; value is the size of the subsets to use")
    parser.add_option("-i", "--skip", default = 1, type = "int", help="interval between snapshots to consider, in whatever units of time snapshots were output in") 
    #parser.add_option("-c", "--correct_formutinf_between_sims", default = "no", type="string", help="correct for excess mutual information between sims")
    #parser.add_option("-l", "--load_matrices_numstructs", default = 0, type = "int", help="if you want to load bootstrap matrices from a previous run, give # of structs per sim")
    #parser.add_option("--plot_2d_histograms", default = False, action = "store_true", help="makes 2d histograms for all pairs of dihedrals in the first bootstrap")
    #parser.add_option("-z", "--zoom_to_step", default = 0, type = "int", help="skips the first n snapshots in xvg files")
    #parser.add_option("-m","--max_num_chis", default = 99, type = "int", help="max number of sidechain chi angles per residue or ligand")
    parser.add_option("-g","--gcc", default = 'gcc', type = "string", help="numpy distutils ccompiler to use. Recommended ones intelem or gcc")
    (options,args)=parser.parse_args()
    mycompiler = options.gcc
    resfile_fn = args[0]
    residue_i = args[1]
    residue_j = args[2]
    prefix = "%s-nsims%d-structs%d-bin%d" % (resfile_fn, options.num_sims, options.num_structs, int(options.binwidth))
    print "prefix: " + str(prefix)
    mymatrix = Mutinf_Matrix_Chis_Bootstraps(prefix)
    Pij_top3, PiPj_top3, top_mutinfs_chi_i, top_mutinfs_chi_j, top_mutinfs, bins, rownames, colnames = mymatrix.pairwise_histogram(residue_i, residue_j, options.xvg_basedir, options.xvg_chidir, options.skip) 
    fig1 = mymatrix.plot_pairwise_histogram((Pij_top3[0])[1:,1:],bins,rownames,colnames)
    fig2 = mymatrix.plot_pairwise_histogram((Pij_top3[1])[1:,1:],bins,rownames,colnames)
    fig3 = mymatrix.plot_pairwise_histogram((Pij_top3[2])[1:,1:],bins,rownames,colnames)
    plt.show()



    
