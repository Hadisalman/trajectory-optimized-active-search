# trajectory-optimized-active-search
Code for the ICRA2018 paper [Trajectory-Optimized Sensing for Active Search of Tissue Abnormalities in Robotic Surgery](https://arxiv.org/pdf/1711.07063.pdf).

The code folder contains two folders; `discrete_probing` and `continuous_probing`.
## discrete_probing
* estimating stiffness map using sequntial discrete probing strategy.
* replicating the results of `Fig.3` in the paper.
#### How to run
* run the `demo.m` file for generating stiffness map using all the algorithms discussed the paper.
* run the `incremental_demo.m` file for generating stiffness map using all the algorithms discussed the paper in an incremental fashion (you can see the points being probed sequentially). This code also saves all the results at after each probe in an automatically generated `results` folder.
* run the `recall_plot_randomized_groundTruths.m` and `recall_plot_same_groundTruth.m` files to generate the recall plot of `Fig. 4`. 


## continuous_probing
* estimating stiffness map using sampling-based trajectory planning to generate continuous palapation paths.
* replicating the results of `Fig.5` in the paper.
#### How to run
* run the `test.m` file for generating stiffness map using any of the algorithms discussed the paper (you can change which algorithm you want to test inside the code) 
* run the `demo.m` file for generating stiffness map using all the algorithms discussed the paper. This code also saves all the results at after each probe in an automatically generated `results` folder.
* run the `recall_plot_randomized_groundTruths.m` and `recall_plot_same_groundTruth.m` files to generate the recall plot of `Fig. 8`. 


# Contact
* Hadi Salman (hadicsalman@gmail.com)
* Elif Ayvali (eayvali@gmail.com)
