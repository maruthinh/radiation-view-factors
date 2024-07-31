# Radiation view factors using the Monte-Carlo Ray-Tracing method

1. This example code uses [PyVista](https://pyvista.org/) package to compute geometric parameters and visualize surfaces, normals, etc.
2. [Numba](https://numba.pydata.org/) JIT compiler is used to improve the performance 
3. An example taken from A Frank, et. al.[^1]. Details of which are shown in the below figure![image](https://github.com/user-attachments/assets/daf07df8-759f-4f27-b680-92c92316ba55)
4. The output of the computed view factors is shown below. We can observe that with 10^6 rays, view factors closely match with the analytical results
```
view factor from surface 1 to 2 is:
  0.04186
time taken is: 1.880399465560913
view factor matrix:
 [[0.       0.043045 0.294282 0.294431 0.183543 0.184128] 
 [0.042769 0.       0.293895 0.295199 0.183512 0.184435] 
 [0.088093 0.087932 0.       0.452062 0.186441 0.186915] 
 [0.087756 0.087865 0.451547 0.       0.185911 0.186542] 
 [0.091866 0.091742 0.312041 0.312048 0.       0.192414] 
 [0.091341 0.092327 0.311122 0.311516 0.19291  0.      ]] 
time taken is:  44.51021385192871
```
5. The geometry considered in the above example is![geom](https://github.com/user-attachments/assets/5dc52d78-4038-4735-bd89-c3ef8ebcc353)
6. The random ray on surface 1 and the normals of surfaces 1 and 2 are shown below![surf12](https://github.com/user-attachments/assets/33f1314a-5149-48a6-9079-e9d6a14c5699)




[^1]: **A Frank, W Heidemann and K Spindler**, Modeling of the surface-to-surface radiation exchange using a Monte Carlo method, *Journal of Physics: Conference Series* 745 (2016) 032143 

