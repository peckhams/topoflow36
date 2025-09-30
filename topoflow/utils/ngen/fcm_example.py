#
#  Copyright, 2025, Scott D. Peckham
#
#  This file is currently in the HARBOR code collection at:
#     github.com/peckhams/topoflow36/topoflow/utils/ngen.
#  Start by cloning the tf36 conda environment as "fuzzy" and then:
#     conda activate fuzzy
#     pip install scikit-fuzzy  # (not avail. for conda install)
#
#  fcm_example.py
#
#------------------------------------------------------------------------
#  To get started:
#
#  % conda activate fuzzy
#  % cd <path-to-fcm_example.py>
#  % python
#  >>> import fcm_example as fcm
#  >>> fcm.run_FCM_example()
#
#------------------------------------------------------------------------
#
#  get_ten_colors()
#
#  run_FCM_example()       # FCM = "Fuzzy C Means"
#  get_FCM_example_data()
#  get_FCM_example_clusters()
#  plot_FCM_example_fpc()
#  get_FCM_example_model()
#  plot_FCM_example_predictions()
#
#------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import skfuzzy as fuzz

#------------------------------------------------------------------------
def get_ten_colors():

    colors = ['b', 'orange', 'g', 'r', 'c', 'm', 'y', 'k',
              'Brown', 'ForestGreen']
    return colors
    
#   get_ten_colors()
#------------------------------------------------------------------------
def run_FCM_example():

    #-------------------------------------------------
    # Note:  Adapted from the example given here:
    # https://scikit-fuzzy.github.io/scikit-fuzzy/
    #    auto_examples/plot_cmeans.html
    # Figures 1 through 5 are generated immediately.
    #-------------------------------------------------
    xpts, ypts, colors = get_FCM_example_data()
    fpcs, alldata = get_FCM_example_clusters( xpts, ypts, colors )
    plot_FCM_example_fpcs( fpcs )

    #-------------------------------------------------    
    # A fuzzy c-means model (FCM) is defined via the
    # cluster centers (and number of clusters)
    #-------------------------------------------------
    cntrs = get_FCM_example_model( alldata )
    
    #-----------------------------------------
    # Generate uniformly sampled data spread
    # across the range [0, 10] in x and y
    #-----------------------------------------
    newdata = np.random.uniform(0, 1, (1100, 2)) * 10
    plot_FCM_example_predictions( newdata, cntrs )
 
#   run_c_means_fuzzy_example()
#------------------------------------------------------------------------
def get_FCM_example_data():

    colors = get_ten_colors()

    #-------------------------------
    # Define three cluster centers
    #-------------------------------
    centers = [[4, 2],
               [1, 7],
               [5, 6]]

    #------------------------------
    # Define three cluster sigmas
    # in x and y, respectively
    #------------------------------
    sigmas = [[0.8, 0.3],
              [0.3, 0.5],
              [1.1, 0.7]]

    #---------------------
    # Generate test data
    #---------------------
    np.random.seed(42)  # Set seed for reproducibility
    xpts = np.zeros(1)
    ypts = np.zeros(1)
    labels = np.zeros(1)
    for i, ((xmu, ymu), (xsigma, ysigma)) in enumerate(zip(centers, sigmas)):
        xpts = np.hstack((xpts, np.random.standard_normal(200) * xsigma + xmu))
        ypts = np.hstack((ypts, np.random.standard_normal(200) * ysigma + ymu))
        labels = np.hstack((labels, np.ones(200) * i))

    #--------------------------
    # Visualize the test data
    #--------------------------
    fig0, ax0 = plt.subplots()
    for label in range(3):
        ax0.plot(xpts[labels == label], ypts[labels == label], '.',
                 color=colors[label])
    ax0.set_title('Test data: 200 points x3 clusters.')

    return xpts, ypts, colors

#   get_FCM_example_data()
#------------------------------------------------------------------------
def get_FCM_example_clusters( xpts, ypts, colors ):

    #---------------------------
    # Set up the loop and plot
    #---------------------------
    fig1, axes1 = plt.subplots(3, 3, figsize=(8, 8))
    alldata = np.vstack((xpts, ypts))
    fpcs = []

    for ncenters, ax in enumerate(axes1.reshape(-1), 2):
        cntrs, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(
            alldata, ncenters, 2, error=0.005, maxiter=1000, init=None)

        #-------------------
        # Store fpc values
        #-------------------
        fpcs.append(fpc)

        #-----------------------------------
        # Plot assigned clusters, for each
        # data point in training set
        #-----------------------------------
        cluster_membership = np.argmax(u, axis=0)
        for j in range(ncenters):
            x = xpts[cluster_membership == j]
            y = ypts[cluster_membership == j]
            ax.plot(x, y, marker='.', linestyle='None', color=colors[j])

        #------------------------------------
        # Mark center of each fuzzy cluster
        #------------------------------------
        for pt in cntrs:
            ax.plot(pt[0], pt[1], 'rs')

        ax.set_title('Centers = {0}; FPC = {1:.2f}'.format(ncenters, fpc))
        ax.axis('off')

    fig1.tight_layout()

    return fpcs, alldata

#   get_FCM_example_clusters()
#------------------------------------------------------------------------
def plot_FCM_example_fpcs( fpcs ):

    fig2, ax2 = plt.subplots()
    ax2.plot(np.r_[2:11], fpcs)
    ax2.set_xlabel("Number of centers")
    ax2.set_ylabel("Fuzzy partition coefficient")
    
#   plot_FCM_example_fpcs()
#------------------------------------------------------------------------
def get_FCM_example_model( alldata, ncenters=3, ndim=2, show=False):

    #-----------------------------------------------------
    # Regenerate fuzzy model with 3 cluster centers.
    # Note: Center ordering is random in this clustering
    # algorithm, so the centers may change places
    #-----------------------------------------------------
    cntrs, u_orig, _, _, _, _, fpc = fuzz.cluster.cmeans(
        alldata, ncenters, 2, error=0.005, maxiter=1000)

#     if (ndim == 3):
#         fig, axes = plt.subplots(3,4, subplot_kw={'projection': '3d'})

    #-----------------------
    # Show 3-cluster model
    #-----------------------
    colors = get_ten_colors()
    for j in range( ncenters ):
        if (ndim == 2):
            fig, ax = plt.subplots()
            x = alldata[0, u_orig.argmax(axis=0) == j]
            y = alldata[1, u_orig.argmax(axis=0) == j]
            ax.plot(x, y, marker='.', linestyle='None',
                    label='series ' + str(j),
                    color=colors[j])
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_xlim(0, 9)
            ax.set_ylim(0, 9)
        else:
            ## fig, ax = plt.subplots( projection='3d')
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            x  = alldata[0, u_orig.argmax(axis=0) == j]
            y  = alldata[1, u_orig.argmax(axis=0) == j]
            z  = alldata[2, u_orig.argmax(axis=0) == j]
            ax.plot(x, y, z, marker='.', linestyle='None',
                    label='series ' + str(j), color=colors[j])
#             axes[j].scatter(alldata[0, u_orig.argmax(axis=0) == j],
#                     alldata[1, u_orig.argmax(axis=0) == j],
#                     alldata[2, u_orig.argmax(axis=0) == j],
#                     marker='.', label='series ' + str(j),
#                     color=colors[j])
            ax.set_xlim(0, 1)
            ax.set_ylim(0,5.6)
            ax.set_zlim(-1,1)
            ax.set_xlabel('Snow fraction')
            ax.set_ylabel('Aridity')
            ax.set_zlabel('Seasonality')
        ax.set_title('Trained model')
        ax.legend()
#         axes[j].set_title('Trained model')
#         axes[j].legend()
     
    print('Number of centers =', ncenters)   
    print('Fuzzy partition coefficient =', fpc)
    if (show):
        plt.show()

    return cntrs
    
#   get_FCM_example_model()
#------------------------------------------------------------------------
def plot_FCM_example_predictions( newdata, cntrs ):

    #-------------------------------------------------------
    # Predict new cluster membership with `cmeans_predict`
    # as well as `cntr` from the 3-cluster model
    #-------------------------------------------------------
    u, u0, d, jm, p, fpc = fuzz.cluster.cmeans_predict(
        newdata.T, cntrs, 2, error=0.005, maxiter=1000)

    #----------------------------------------------
    # Harden cluster membership for visualization
    #----------------------------------------------
    cluster_membership = np.argmax(u, axis=0)

    #---------------------------------------------------------------
    # Plot the classified uniform data.  Note for visualization
    # the maximum membership value has been taken at each point
    # (i.e. these are hardened, not fuzzy results visualized)
    # but the full fuzzy result is the output from cmeans_predict.
    #---------------------------------------------------------------
    fig3, ax3 = plt.subplots()
    ax3.set_title('Random points classified according to known centers')
    for j in range(3):
        ax3.plot(newdata[cluster_membership == j, 0],
                 newdata[cluster_membership == j, 1], 'o',
                 label='series ' + str(j))
    ax3.legend()
    plt.show()

#   plot_FCM_example_predictions()
#------------------------------------------------------------------------
