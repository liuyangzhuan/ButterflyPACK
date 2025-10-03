make
export GMSH_DIR=/global/homes/l/liuyangz/Cori/my_software/gmsh-4.6.0-Linux64/bin
#export GMSH_DIR=/home/administrator/Desktop/research/preprocessor/gmsh-2.11.0-Linux/bin
# export GMSH_DIR=/home/administrator/Desktop/software/gmsh-2.11.0-Linux/bin
# export GMSH_DIR=/home/administrator/Desktop/Software/gmsh-4.6.0-Linux64/bin
# export GMSH_DIR=/Applications/Gmsh.app/Contents/MacOS/gmsh/bin
# # # # #${GMSH_DIR}/gmsh sphere.geo -1 -2 -format msh -algo del2d -clmin 7.65e-1 -clmax 1e0   # 48 patches
# # # # #${GMSH_DIR}/gmsh sphere.geo -1 -2 -format msh -algo del2d -clmin 3.8e-1 -clmax 5e-1   # 192 patches
# # # # ##${GMSH_DIR}/gmsh sphere.geo -1 -2 -format msh -algo del2d -clmin 1.9e-1 -clmax 2.5e-1   # 747 patches
# # # # # ${GMSH_DIR}/gmsh sphere.geo -1 -2 -format msh -algo del2d -clmin 0.95e-1 -clmax 1.25e-1   # 2427 patches
# # # # #${GMSH_DIR}/gmsh sphere.geo -1 -2 -format msh -algo del2d -clmin 0.475e-1 -clmax 0.625e-1   # 9487 patches
# # # # # ${GMSH_DIR}/gmsh sphere.geo -2 -format msh -algo del2d -clmin 1.2e-2 -clmax 1.56e-2   # 137352 patches
# # # # # ${GMSH_DIR}/gmsh sphere.geo -2 -format msh -algo del2d -clmin 0.6e-2 -clmax 0.78e-2   # 552655 patches
# # # # # ${GMSH_DIR}/gmsh sphere.geo -2 -format msh -algo del2d -clmin 0.3e-2 -clmax 0.39e-2 #2189137 patches
# # # # # ${GMSH_DIR}/gmsh sphere.geo -2 -format msh -algo del2d -clmin 0.15e-2 -clmax 0.2e-2 #8331149 patches
##


# name=corner_2500
# ${GMSH_DIR}/gmsh corner.geo -2 -o $name.nas -algo del2d -clmin 0.25e-1 -clmax 0.4e-1 -string "Mesh.BdfFieldFormat = 2;"  
# ./con_all ${name}.nas
# mv node.geo ${name}_node.inp
# mv elem.geo ${name}_elem.inp

# name=corner_10000
# ${GMSH_DIR}/gmsh corner.geo -2 -o $name.nas -algo del2d -clmin 0.125e-1 -clmax 0.2e-1 -string "Mesh.BdfFieldFormat = 2;"  
# ./con_all ${name}.nas
# mv node.geo ${name}_node.inp
# mv elem.geo ${name}_elem.inp

# name=corner_40000
# ${GMSH_DIR}/gmsh corner.geo -2 -o $name.nas -algo del2d -clmin 0.0625e-1 -clmax 0.1e-1 -string "Mesh.BdfFieldFormat = 2;"  
# ./con_all ${name}.nas
# mv node.geo ${name}_node.inp
# mv elem.geo ${name}_elem.inp



# name=halfsphere_1200
# ${GMSH_DIR}/gmsh halfsphere.geo -2 -o $name.nas -algo del2d -clmin 0.8e-1 -clmax 1.3e-1 -string "Mesh.BdfFieldFormat = 2;"  
# ./con_all ${name}.nas
# mv node.geo ${name}_node.inp
# mv elem.geo ${name}_elem.inp


# name=halfsphere_2300
# ${GMSH_DIR}/gmsh halfsphere.geo -2 -o $name.nas -algo del2d -clmin 0.6e-1 -clmax 0.9e-1 -string "Mesh.BdfFieldFormat = 2;"  
# ./con_all ${name}.nas
# mv node.geo ${name}_node.inp
# mv elem.geo ${name}_elem.inp

# name=halfsphere_9000
# ${GMSH_DIR}/gmsh halfsphere.geo -2 -o $name.nas -algo del2d -clmin 0.3e-1 -clmax 0.44e-1 -string "Mesh.BdfFieldFormat = 2;"  
# ./con_all ${name}.nas
# mv node.geo ${name}_node.inp
# mv elem.geo ${name}_elem.inp

# name=halfsphere_32000
# ${GMSH_DIR}/gmsh halfsphere.geo -2 -o $name.nas -algo del2d -clmin 0.15e-1 -clmax 0.22e-1 -string "Mesh.BdfFieldFormat = 2;"  
# ./con_all ${name}.nas
# mv node.geo ${name}_node.inp
# mv elem.geo ${name}_elem.inp


# ${GMSH_DIR}/gmsh halfsphere.geo -2 -o halfsphere_35000.nas -algo del2d -clmin 0.15e-1 -clmax 0.22e-1 -string "Mesh.BdfFieldFormat = 2;"
#

# name=sphere_1200
# ${GMSH_DIR}/gmsh sphere.geo -2 -o $name.nas -algo del2d -clmin 1.2e-1 -clmax 1.8e-1 -string "Mesh.BdfFieldFormat = 2;"  
# ./con_all ${name}.nas
# mv node.geo ${name}_node.inp
# mv elem.geo ${name}_elem.inp


#  name=sphere_2300
#  ${GMSH_DIR}/gmsh sphere.geo -2 -o $name.nas -algo del2d -clmin 0.8e-1 -clmax 1.3e-1 -string "Mesh.BdfFieldFormat = 2;"  
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp

#  name=sphere_9000
#  ${GMSH_DIR}/gmsh sphere.geo -2 -o $name.nas -algo del2d -clmin 0.4e-1 -clmax 0.65e-1 -string "Mesh.BdfFieldFormat = 2;"  
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp

#  name=sphere_32000
#  ${GMSH_DIR}/gmsh sphere.geo -2 -o $name.nas -algo del2d -clmin 0.2e-1 -clmax 0.325e-1 -string "Mesh.BdfFieldFormat = 2;"  
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp

#  name=sphere_128000
#  ${GMSH_DIR}/gmsh sphere.geo -2 -o $name.nas -algo del2d -clmin 0.1e-1 -clmax 0.1625e-1 -string "Mesh.BdfFieldFormat = 2;"  
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp


#  name=sphere_512000
#  ${GMSH_DIR}/gmsh sphere.geo -2 -o $name.nas -algo del2d -clmin 0.05e-1 -clmax 0.08125e-1 -string "Mesh.BdfFieldFormat = 2;"  
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp


 name=sphere_2M
#  ${GMSH_DIR}/gmsh sphere.geo -2 -o $name.nas -algo del2d -clmin 0.025e-1 -clmax 0.040625e-1 -string "Mesh.BdfFieldFormat = 2;"  
 ./con_all ${name}.nas
 mv node.geo ${name}_node.inp
 mv elem.geo ${name}_elem.inp


 #name=sphere_17999
 #${GMSH_DIR}/gmsh sphere.geo -2 -o $name.nas -algo del2d -clmin 0.3e-1 -clmax 0.44e-1 -string "Mesh.BdfFieldFormat = 2;"  
 #./con_all ${name}.nas
 #mv node.geo ${name}_node.inp
 #mv elem.geo ${name}_elem.inp

 #name=sphere_71000
 #${GMSH_DIR}/gmsh sphere.geo -2 -o $name.nas -algo del2d -clmin 0.15e-1 -clmax 0.22e-1 -string "Mesh.BdfFieldFormat = 2;"  
 #./con_all ${name}.nas
 #mv node.geo ${name}_node.inp
 #mv elem.geo ${name}_elem.inp


#${GMSH_DIR}/gmsh sphere.geo -2 -o sphere_4500.nas -algo del2d -clmin 0.6e-1 -clmax 0.9e-1 -string "Mesh.BdfFieldFormat = 2;"  
#${GMSH_DIR}/gmsh sphere.geo -2 -o sphere_17999.nas -algo del2d -clmin 0.3e-1 -clmax 0.44e-1 -string "Mesh.BdfFieldFormat = 2;"
#${GMSH_DIR}/gmsh sphere.geo -2 -o sphere_71000.nas -algo del2d -clmin 0.15e-1 -clmax 0.22e-1 -string "Mesh.BdfFieldFormat = 2;"
#
# name=plate_1000
# ${GMSH_DIR}/gmsh plate.geo -2 -o $name.nas -algo front2d -clmin 1.0e-1 -clmax 1.0e-1  -string "Mesh.BdfFieldFormat = 2;" # 583016 patches
# ./con_all ${name}.nas
# mv node.geo ${name}_node.inp
# mv elem.geo ${name}_elem.inp

# name=plate_2000
# ${GMSH_DIR}/gmsh plate.geo -2 -o $name.nas -algo front2d -clmin 0.73e-1 -clmax 0.73e-1  -string "Mesh.BdfFieldFormat = 2;" # 1946 patches
# ./con_all ${name}.nas
# mv node.geo ${name}_node.inp
# mv elem.geo ${name}_elem.inp

# name=plate_8000
# ${GMSH_DIR}/gmsh plate.geo -2 -o $name.nas -algo front2d -clmin 0.36e-1 -clmax 0.36e-1  -string "Mesh.BdfFieldFormat = 2;" # 7552 patches
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp

#  name=plate_32000
#  ${GMSH_DIR}/gmsh plate.geo -2 -o $name.nas -algo front2d -clmin 0.18e-1 -clmax 0.18e-1  -string "Mesh.BdfFieldFormat = 2;" # 7552 patches
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp

#  name=plate_128000
#  ${GMSH_DIR}/gmsh plate.geo -2 -o $name.nas -algo front2d -clmin 0.09e-1 -clmax 0.09e-1  -string "Mesh.BdfFieldFormat = 2;" # 7552 patches
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp

#  name=plate_512000
#  ${GMSH_DIR}/gmsh plate.geo -2 -o $name.nas -algo front2d -clmin 0.045e-1 -clmax 0.045e-1  -string "Mesh.BdfFieldFormat = 2;" # 7552 patches
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp


#name=plate_10240000
#${GMSH_DIR}/gmsh plate.geo -2 -o $name.nas -algo front2d -clmin 0.01125e-1 -clmax 0.01125e-1  -string "Mesh.BdfFieldFormat = 2;" # 7552 patches
#./con_all ${name}.nas
#mv node.geo ${name}_node.inp
#mv elem.geo ${name}_elem.inp


#${GMSH_DIR}/gmsh plate.geo -2 -o plate_580000.nas -algo front2d -clmin 0.4e-2 -clmax 0.4e-2   # 583016 patches
# ${GMSH_DIR}/gmsh plate.geo -2 -o plate_2300000.nas -algo front2d -clmin 0.2e-2 -clmax 0.2e-2   # 2314782 patches
#${GMSH_DIR}/gmsh plate.geo -2 -o plate_10000000.nas -algo front2d -clmin 0.1e-2 -clmax 0.1e-2   # 2314782 patches


#  name=pillbox_1000
#  ${GMSH_DIR}/gmsh pillbox.geo -2 -o $name.nas -algo del2d -clmin 0.12e-1 -clmax 0.2e-1 -string "Mesh.BdfFieldFormat = 2;"  
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp


#  name=pillbox_4000
#  ${GMSH_DIR}/gmsh pillbox.geo -2 -o $name.nas -algo del2d -clmin 0.06e-1 -clmax 0.1e-1 -string "Mesh.BdfFieldFormat = 2;"  
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp



#  name=pillbox_50K
#  ${GMSH_DIR}/gmsh pillbox.geo -2 -o $name.nas -algo del2d -clmin 0.015e-1 -clmax 0.025e-1 -string "Mesh.BdfFieldFormat = 2;"  
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp


#  name=pillbox_200K
#  ${GMSH_DIR}/gmsh pillbox.geo -2 -o $name.nas -algo del2d -clmin 0.0075e-1 -clmax 0.0125e-1 -string "Mesh.BdfFieldFormat = 2;"  
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp


#  name=cavity_wakefield_4K_feko
# #  ${GMSH_DIR}/gmsh cavity_wakfield.geo -2 -o $name.nas -algo del2d -clmin 0.12e2 -clmax 0.2e2 -string "Mesh.BdfFieldFormat = 2;"  
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp

#  name=cavity_wakefield_14K_feko
# #  ${GMSH_DIR}/gmsh cavity_wakfield.geo -2 -o $name.nas -algo del2d -clmin 0.12e2 -clmax 0.2e2 -string "Mesh.BdfFieldFormat = 2;"  
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp

#  name=cavity_rec_5K_feko
# #  ${GMSH_DIR}/gmsh cavity_wakfield.geo -2 -o $name.nas -algo del2d -clmin 0.12e2 -clmax 0.2e2 -string "Mesh.BdfFieldFormat = 2;"  
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp

#  name=cavity_rec_17K_feko
# #  ${GMSH_DIR}/gmsh cavity_wakfield.geo -2 -o $name.nas -algo del2d -clmin 0.12e2 -clmax 0.2e2 -string "Mesh.BdfFieldFormat = 2;"  
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp


#  name=cavity_5cell_400K_feko
# #  ${GMSH_DIR}/gmsh cavity_wakfield.geo -2 -o $name.nas -algo del2d -clmin 0.12e2 -clmax 0.2e2 -string "Mesh.BdfFieldFormat = 2;"  
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp


#  name=cavity_5cell_1M_feko
# #  ${GMSH_DIR}/gmsh cavity_wakfield.geo -2 -o $name.nas -algo del2d -clmin 0.12e2 -clmax 0.2e2 -string "Mesh.BdfFieldFormat = 2;"  
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp


#  name=cavity_5cell_400K_feko
# #  ${GMSH_DIR}/gmsh cavity_wakfield.geo -2 -o $name.nas -algo del2d -clmin 0.12e2 -clmax 0.2e2 -string "Mesh.BdfFieldFormat = 2;"  
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp

#  name=cavity_5cell_30K_feko
# #  ${GMSH_DIR}/gmsh cavity_wakfield.geo -2 -o $name.nas -algo del2d -clmin 0.12e2 -clmax 0.2e2 -string "Mesh.BdfFieldFormat = 2;"  
#  ./con_all ${name}.nas
#  mv node.geo ${name}_node.inp
#  mv elem.geo ${name}_elem.inp


#  rm *.nas

#
