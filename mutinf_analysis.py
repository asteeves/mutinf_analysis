#!/sw/bin/python


import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
from mpl_toolkits.mplot3d import Axes3D
from numpy import float64
import AnnoteModule

class mutInfmat(object):
    """
    A class for holding the information associated with a mutual information matrix.
    
    In typical usage, once the  mutual information matrix is instantiated,
    a group of 2D plots of the eigenvectors is examined (twoDplots),
    the best 2D plot is used to 
    
    """
    def __init__(self, myfilename):
        self.mymatrix, self.resnames, self.junk  = read_res_matrix(myfilename)
        self.eigvals, self.eigvecs = sortedeig(self.mymatrix)
    
    def unsort(self):
        """
	Plots the unsorted mutual information matrix.
	"""
        fig = plt.figure(facecolor='w')
	ax = fig.add_subplot(111)
	im = ax.matshow(self.mymatrix,cmap=plt.cm.Blues)
	im.set_clim(0,0.5)
	ax.set_yticks(np.arange(0,len(self.resnames)))
	ax.set_yticklabels(self.resnames)
	ax.set_xticks([])
	for tick in ax.yaxis.get_major_ticks():
	    tick.label1On = False
	    tick.label2On = True
	    tick.tick1On = False
	    tick.tick2On = False
	for label in ax.yaxis.get_ticklabels():
	    label.set_size(4)
	ax.set_title('Unclustered MutInf matrix')
	fig.show()
	return fig

    def heatmap(self):
        """
	Plots the mutual information matrix as a heatmap with the associated dendrogram.
	Hierarchical clustering uses the default parameters of the hclust function in R.
	(i.e. 'complete' clustering with a euclidean distance metric)
	"""
        p = ssd.pdist(self.mymatrix)
        links = sch.linkage(p,method='complete')
        dend = sch.dendrogram(links,no_plot=True)
        idx = dend['leaves']
	#reorder the input matrices and labels according to the index determined by the clustering
        reordered = self.mymatrix[idx,:]
        reordered = reordered[:,idx]
	reordres = [self.resnames[i] for i in idx]
        fig = plt.figure(figsize=(8,10),facecolor='w')
        ax1 = fig.add_axes([0.1,0.85,0.8,0.08],frameon=False)
        dend = sch.dendrogram(links,orientation='top')
        ax1.set_xticks([])
        ax1.set_yticks([])
        axmatrix = fig.add_axes([0.1,0.1,0.8,0.7])
        idx = dend['leaves']
        reordered = self.mymatrix[idx,:]
        reordered = reordered[:,idx]
        im = axmatrix.matshow(reordered, aspect='auto', origin='upper', cmap=plt.cm.binary)
        im.set_clim(0,1)
        axmatrix.set_xticks(np.arange(1,len(self.resnames)))
        axmatrix.set_xticklabels(reordres,rotation='vertical')
        axmatrix.set_yticks([])
        for tick in axmatrix.xaxis.get_major_ticks():
	    tick.tick1On = False
	    tick.tick2On = False
	for label in axmatrix.xaxis.get_ticklabels():
	    label.set_rotation(90)
	    label.set_size(6)
	fig.suptitle('Clustered MutInf Heatmap',size=12.0)
        #fig.show()

    def twoDplots(self):
        """
        Provides a plot for a overview of the two-dimensional slices of eigenvector space.
	TODO:Should be configurable for which 2D plots are desired.
	"""
        fig = plt.figure(figsize=(12,4),facecolor='w')
	s1 = fig.add_subplot(131)
	s1.scatter(self.eigvecs[:,0],self.eigvecs[:,1])
	s1.set_xlabel(r'$|0>$')
	s1.set_ylabel(r'$|1>$')
	s2 = fig.add_subplot(132)
	s2.scatter(self.eigvecs[:,1],self.eigvecs[:,2])
	s2.set_xlabel(r'$|1>$')
	s2.set_ylabel(r'$|2>$')
	s3 = fig.add_subplot(133)
	s3.scatter(self.eigvecs[:,1],self.eigvecs[:,3])
	s3.set_xlabel(r'$|1>$')
	s3.set_ylabel(r'$|3>$')
	fig.suptitle('Eigenvector projections')
        return fig
#	fig.show()

    def threeDplot(self,ind1,ind2,ind3):
        fig = plt.figure()
	ax = Axes3D(fig)
        ax.scatter(self.eigvecs[:,ind1],self.eigvecs[:,ind2],self.eigvecs[:,ind3])
        return fig
        #fig.show()
    
    def twoDxcc(self,eigind1,eigind2):
        """
	Scan a trial ray around in a 2D space 
	"""	
	mymaps = (plt.cm.Blues, plt.cm.Reds, plt.cm.Greens)
	eigv1=self.eigvecs[:,eigind1]
	eigv2=self.eigvecs[:,eigind2]
	numres=self.mymatrix.shape[0]
	steps = 360
	alpha= np.zeros(steps)
	meritfunc = np.zeros(steps)
	poscontrib = np.zeros((steps,numres))

	#sigma0 currently a standin for a noise determined from the data, determines how many 
	#sectors are identified in 
	sigma0 = 0.1
	#similar cutoff determines the magnitude of the contribution required for sector membership
	cutoff = 0.1 
	
	for ai in np.arange(steps):
	    alpha[ai] = ai*(360/steps)*np.pi/180
            #print angle, radangle, np.sin(radangle), np.cos(radangle)
	    for resi in np.arange(numres):
	        #calculate value of the merit function summed over all residues
                proj = eigv1[resi]*np.cos(alpha[ai])+eigv2[resi]*np.sin(alpha[ai])
		#dist = eigv2[resi]/np.cos(alpha[ai]) - eigv1[resi]*np.sin(alpha[ai])/(np.cos(alpha[ai]))**2
		dist = np.sqrt(eigv1[resi]**2+eigv2[resi]**2-proj**2)
		contrib =proj*np.exp(-0.5*(dist/sigma0)**2)
		if contrib>0:
		    poscontrib[ai,resi]=contrib
		    meritfunc[ai]=meritfunc[ai]+contrib
	
	maxlocs, maxvals  = extrema(meritfunc,min=False)
	maxangles = alpha[maxlocs]

	fig = plt.figure(figsize=(8,4))
	ax1=fig.add_subplot(121)
        ax1.scatter(eigv1,eigv2,c='w',marker='o')
       	ax2 = fig.add_subplot(122)
	ax2.plot(alpha*360/(2*np.pi),meritfunc,'k-')
	ax2.stem(maxangles*360/(2*np.pi),maxvals,linefmt='r-',markerfmt='k+')
	ax2.set_xlim(0,360);ax2.set_xlabel('Degrees')
	ax2.set_ylabel('Merit Function')
	ax2.set_xticks([0,90,180,270,360])
	for ray in maxangles:
	    vecx = [0, np.cos(ray)]
	    vecy = [0, np.sin(ray)]
	    ax1.plot(vecx,vecy,'r-')
        ax1.set_xlim(-0.5,0.5);ax1.set_xlabel('|'+str(eigind1)+'>')
	ax1.set_ylim(-0.5,0.5);ax1.set_ylabel('|'+str(eigind2)+'>')
        fig.suptitle('Sector analysis of MutInfMat using eigs '+str(eigind1)+' and '+str(eigind2),
	              size=12.0)

	listofsecs = []

	for secnum, secpos in enumerate(maxlocs):
	    #find the residues that survive the cutoff criteria for each sector
	    secname = 'sector'+str(secnum)
	    secind = np.nonzero(poscontrib[secpos,:]>cutoff)[0]#to get array from tuple
	    weights = poscontrib[secpos,secind]
	    seccol = mymaps[secnum]

            current = sector(secname,secind,weights,seccol)
            listofsecs.append(current)

	    #ax1.scatter(eigv1[secind],eigv2[secind],c=poscontrib[secpos,secind],cmap=mymaps[secnum])

	#return maxlocs, maxvals, maxangles, poscontrib
	self.listofsecs=listofsecs

    def twoDplot_vec(self,limit):
        """
	Provides an array of plots of 2-D slices in the space of the eigenvectors of the 
	mutual information matrix. Points are colored according to sectors (if any), 
	with the sector contributions determining the saturation.
	
	These plots have linked interactive annotations, requiring the AnnotationModule.
	"""
        fig=plt.figure(figsize=(12,12),facecolor='w')
	grid=limit-1
	afs=[]
        for vec1 in range(limit):
	    for vec2 in np.arange(vec1+1,limit):
	        index = grid*(vec1)+(vec2)
                s1=fig.add_subplot(grid,grid,index)
		s1.scatter(self.eigvecs[:,vec1],self.eigvecs[:,vec2],c='w')
	        #s1.xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=3))	
	        #s1.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=3))	
		afs.append(AnnoteModule.AnnoteFinder(self.eigvecs[:,vec1],
		                 self.eigvecs[:,vec2],self.resnames,xtol=5,ytol=5))
		plt.connect('button_press_event',afs[-1])
		for sec in self.listofsecs:
		    s1.scatter(self.eigvecs[sec.members,vec1],self.eigvecs[sec.members,vec2],
		              c=sec.weights,cmap=sec.seccolor,
			      norm=matplotlib.colors.Normalize(vmin=0,vmax=0.3))    
        #using linked annotation finders from MatplotLib/Interactive Plotting
        for i in range(len(afs)):
	    allButSelfAfs = afs[:i]+afs[i+1:]
	    afs[i].links.extend(allButSelfAfs)


    def filteredmat(self,eiglist):
        """
        Filters the mutual information matrix, keeping only the eigenvectors given in the list.
        """
        #reconstruct matrix from the eigenvectors given in the eiglist, (sca spectral cleaning)
	filtmat = np.zeros(self.mymatrix.shape)
        for eigind in eiglist:
	    #not well tested, use of np.outer
	    filtmat = filtmat+self.eigvals[eigind]*np.outer(self.eigvecs[:,eigind],self.eigvecs[:,eigind])
	return filtmat


def extrema(x, max = True, min = True, strict = False, withend = False):
    """
    This function will index the extrema of a given array x.
	
    Options:
		max		If true, will index maxima
		min		If true, will index minima
		strict		If true, will not index changes to zero gradient
		withend	If true, always include x[0] and x[-1]
	
    This function will return a tuple of extrema indicies and values
    """
	
    # This is the gradient
    from numpy import zeros
    dx = zeros(len(x))
    from numpy import diff
    dx[1:] = diff(x)
    dx[0] = dx[1]
	
    # Clean up the gradient in order to pick out any change of sign
    from numpy import sign
    dx = sign(dx)
	
    # define the threshold for whether to pick out changes to zero gradient
    threshold = 0
    if strict:
    	threshold = 1
		
    # Second order diff to pick out the spikes
    d2x = diff(dx)
	
    if max and min:
     	d2x = abs(d2x)
    elif max:
 	d2x = -d2x
	
    # Take care of the two ends
    if withend:
	d2x[0] = 2
 	d2x[-1] = 2
	
    # Sift out the list of extremas
    from numpy import nonzero
    ind = nonzero(d2x > threshold)[0]
	
    return ind, x[ind]

class sector(object):
    """
    a sector can be defined as in SCA 4.0, or as a direction in the space spanned by 
    the eigenvectors associated with the largest (most significant) eigenvalues
    """
    def __init__(self,name,members,weights,seccolor):
        self.name = name
	self.members = members
        self.weights = weights
        self.seccolor = seccolor

def sortedeig(mymatrix):
    """
    Diagonalizes the matrix given to it, and returns the eigenvectors and eigenvalues
    sorted according to decreasing magnitude of the eigenvalue (largest eigenvalue first).
    """
    eigvals, eigvecs = np.linalg.eigh(mymatrix)
    indx = np.argsort(np.abs(eigvals))
    #resort eigenvalues and eigenvectors
    eigvals = eigvals[indx[::-1]] # reverses the order of the index 
    eigvecs = eigvecs[:,indx[::-1]]
    return eigvals, eigvecs

def read_res_matrix(myfilename):
    """
    very slightly modified version of read_res_matrix from res_utils.py.
    Changed appended name in myname_num loop.
    """
    rownames = []
    colnames = []
    myfile = open(myfilename, 'r')
    inlines = myfile.readlines()
    myfile.close()
    res = inlines[0].split()
    mymatrix = np.zeros((int(len(inlines[1:])), int(len(res))), float64)

    for myname_num in res:
        colnames.append(myname_num)
    for row_num in range(int(len(inlines[1:]))):
        thisline = inlines[row_num + 1]
        thislinedata = thisline.split()
        thisname = thislinedata[0]
        res_num = int(np.floor(row_num))
        thislinenums = map(float, thislinedata[1:])
        thislinearray = np.array(thislinenums, float64)
        rownames.append(thisname)
        for col_num in range(len(colnames)):
            mymatrix[res_num,col_num] = float64(thislinearray[col_num])
    return mymatrix, rownames, colnames

if __name__ == "__main__":
    print 'You are using matplotlib version '+matplotlib.__version__
    print 'You can make changes to global plotting options in '+matplotlib.matplotlib_fname()
    j = mutInfmat('2esk_demo.txt')
    fig1 = j.twoDxcc(1,2)
    fig2 = j.twoDplot_vec(4)
    fig3 = j.threeDplot(0,1,2)
    plt.show()
