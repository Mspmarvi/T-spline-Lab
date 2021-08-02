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


This code is to build a T-spline surface by MATLAB

how to build a T-spline surface:

step 1: build a 3 dimensional matrix with size (n+1)*(m+1)*4  control mesh named (cp if you want), the third dimension means (X,Y,Z,weight) coordinate

step 2: Name_of_surface = T_Surface(cp), and now you build a T-spline NURBS surface in T-spine format

there is an explanation in code for each functions

examples:

![image](https://user-images.githubusercontent.com/34415658/124386176-04dd2280-dc8e-11eb-8f04-52e2fa764b4f.png)
![image](https://user-images.githubusercontent.com/34415658/127885732-66ff8627-d2d8-4692-a577-c0988306cbd8.png)

![image](https://user-images.githubusercontent.com/34415658/124386179-09094000-dc8e-11eb-961c-32e9ee69933d.png)
![image](https://user-images.githubusercontent.com/34415658/127885759-8bf797a8-2b87-4e24-9035-3e15312529af.png)

![image](https://user-images.githubusercontent.com/34415658/124386188-0c9cc700-dc8e-11eb-9c30-466714480c55.png)
![image](https://user-images.githubusercontent.com/34415658/127886692-f226b90f-b8b5-45bf-9a04-2b13017e3ab2.png)

![image](https://user-images.githubusercontent.com/34415658/124386197-10304e00-dc8e-11eb-9d32-45356754ad75.png)
![image](https://user-images.githubusercontent.com/34415658/127886729-992c99d7-58b7-439b-b75f-ec746f47ed44.png)

![image](https://user-images.githubusercontent.com/34415658/127886966-043eff05-c2e0-4838-9964-e07cbdadb0bb.png)
![image](https://user-images.githubusercontent.com/34415658/127887007-33008c1a-eaa1-42a2-9ef1-e2f3251df75d.png)
![image](https://user-images.githubusercontent.com/34415658/127887024-623d6880-4474-412a-b446-e2a53bb47f2f.png)
![image](https://user-images.githubusercontent.com/34415658/127887033-6bff3ecb-ab27-4143-bc4a-ce9dab6fead9.png)

![image](https://user-images.githubusercontent.com/34415658/124390081-e9c6de80-dc9e-11eb-930b-1b5e2c898bac.png)

