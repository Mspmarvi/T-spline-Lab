# T-spline Lab                                 created by Mohammad Sadeghi Pour (Email: Mspmarvi@gmail.com)
T-Spline Lab is an object oriented package to build and manipulate T-spline surfaces

abilities 

1- building T-spline surface

2- inserting contol points based on approach mentioned in "T-splines and T-NURCCs" 2003

3- inserting contol points based on local refinement algorithm

4- merging T-spline surfaces one by one with continuity C0, C1 and C2

5- plotting the topology of control mesh

6- plotting the parametric view of control mesh

7- inserting a row or collumn of control points

8- you can change manually every single objects of T-spline surface including control point properties, edge properties, uknot and vknot properties

9- and there are a lots of functions to help you out for your researches


This code is to build a T-spline surfac eby MATLAB

how to build a T-spline surface:

step 1: build a 3 dimensional matrix with size (n+1)*(m+1)*4  control mesh named (cp if you want), the third dimension means (X,Y,Z,weight) coordinate

step 2: Name_of_surface = T_Surface(cp), and now you build a T-spline NURBS surface in T-spine format

there is an explanation in code for each functions

examples:

![image](https://user-images.githubusercontent.com/34415658/124386176-04dd2280-dc8e-11eb-8f04-52e2fa764b4f.png)![image](https://user-images.githubusercontent.com/34415658/124386179-09094000-dc8e-11eb-961c-32e9ee69933d.png)
![image](https://user-images.githubusercontent.com/34415658/124386188-0c9cc700-dc8e-11eb-9c30-466714480c55.png)![image](https://user-images.githubusercontent.com/34415658/124386197-10304e00-dc8e-11eb-9d32-45356754ad75.png)
![image](https://user-images.githubusercontent.com/34415658/124390081-e9c6de80-dc9e-11eb-930b-1b5e2c898bac.png)


