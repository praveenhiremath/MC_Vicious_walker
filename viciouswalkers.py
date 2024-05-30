## Author: Praveenkumar Hiremath (Lund University 10/20/2019.)
from __future__ import division
import numpy as np
import matplotlib
from matplotlib import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math 
from math import *
import io
import os
from pylab import *
from scipy.optimize import leastsq
from matplotlib.pyplot import figure
figure(num=None, figsize=(8, 20), dpi=80, facecolor='w', edgecolor='k')
rcParams['legend.fontsize'] = 10
import random
from random import *

'''
if len(sys.argv)==2:
  out_file = sys.argv[1]
else:
  print "Usage: python", sys.argv[0], "output_file_name"
  sys.exit(0)
'''

def get_inputs():
   Dim=input("Enter the dimension of the problem\n1-for 1-D, 2-for 2-D, 3- for 3D\n")
   if (Dim==1):
     print "1D case: Probability of move in all directions (right/left) is equal."
     probabilities=np.zeros(shape=(2),dtype=float)     
     num_walkers=input("Enter the number of walkers\n")
     probabilities[0]=0.5
     probabilities[1]=0.5

   if (Dim==2):
     print "2D case: Probability of move in all directions (up/down/left/right) is equal."
     probabilities=np.zeros(shape=(4),dtype=float)
     num_walkers=input("Enter the number of walkers\n")
     move_prob_x=0.25
     move_prob_nx=move_prob_x
     move_prob_y=move_prob_x
     move_prob_ny=move_prob_x
     probabilities[0]=move_prob_x
     probabilities[1]=move_prob_nx
     probabilities[2]=move_prob_y
     probabilities[3]=move_prob_ny

   if (Dim==3):
     print "3D case: Probability of move in all directions (+ve and -ve directions along x,y,z axes) is equal."
     probabilities=np.zeros(shape=(6),dtype=float)
     num_walkers=input("Enter the number of walkers\n")
     move_prob_x=1/6
     move_prob_nx=move_prob_x
     move_prob_y=move_prob_x
     move_prob_ny=move_prob_x
     move_prob_z=move_prob_x
     move_prob_nz=move_prob_x
     probabilities[0]=move_prob_x
     probabilities[1]=move_prob_nx
     probabilities[2]=move_prob_y
     probabilities[3]=move_prob_ny
     probabilities[4]=move_prob_z
     probabilities[5]=move_prob_nz

   return(int(Dim),int(num_walkers), probabilities)

    
def sim_grid(Dim):
  if (Dim==1):
     x_length=input("Enter length in x-direction\n")
     num_pts_x=input("Enter number of grid points along x\n")
     num_pts_y=0
     num_pts_z=0
     xyz_grid_pts=np.zeros(shape=(num_pts_x),dtype=int)
     rand_i=np.linspace(0,num_pts_x-10,num_walkers).astype(int)
     np.random.shuffle(rand_i)
     rand_j=np.zeros(shape=(num_walkers),dtype=int)
     rand_k=np.zeros(shape=(num_walkers),dtype=int)
     

  if (Dim==2):
     x_length=input("Enter length in x-direction\n")
     y_length=input("Enter length in y-direction\n")
     num_pts_x=input("Enter number of grid points along x\n")
     num_pts_y=input("Enter number of grid points along y\n")
     num_pts_z=0
     xyz_grid_pts=np.zeros(shape=(num_pts_x,num_pts_y),dtype=int)
     rand_i=np.linspace(0,num_pts_x-10,num_walkers).astype(int)
     np.random.shuffle(rand_i)
     rand_j=np.linspace(0,num_pts_y-10,num_walkers).astype(int)
     np.random.shuffle(rand_j)
     rand_k=np.zeros(shape=(num_walkers),dtype=int)

  if (Dim==3):
     x_length=input("Enter length in x-direction\n")
     y_length=input("Enter length in y-direction\n")
     z_length=input("Enter length in z-direction\n")
     num_pts_x=input("Enter number of grid points along x\n")
     num_pts_y=input("Enter number of grid points along y\n")
     num_pts_z=input("Enter number of grid points along z\n")
     xyz_grid_pts=np.zeros(shape=(num_pts_x,num_pts_y,num_pts_z),dtype=int)
     rand_i=np.linspace(0,num_pts_x-10,num_walkers).astype(int)
     np.random.shuffle(rand_i)
     rand_j=np.linspace(0,num_pts_y-10,num_walkers).astype(int)
     np.random.shuffle(rand_j)
     rand_k=np.linspace(0,num_pts_z-10,num_walkers).astype(int)
     np.random.shuffle(rand_k)

  return rand_i, rand_j, rand_k, xyz_grid_pts,num_pts_x,num_pts_y,num_pts_z

Dim, num_walkers, probabilities = get_inputs()
print "Problem is ", int(Dim), "number of walkers/particles = ",int(num_walkers),"with equal probability to move to neighbouring sites\n"

rand_i, rand_j, rand_k, xyz_grid_pts,num_pts_x,num_pts_y,num_pts_z=sim_grid(Dim)
#print "random initial positions of walkers", rand_i, "\n"
print "size of xyz_grid_pts",xyz_grid_pts.shape

def create_random_particles(Dim,xyz_grid_pts,num_walkers,rand_i,rand_j,rand_k,probabilities):
  if (Dim==1):
     print "number of walkers",num_walkers
     for j in range(0,num_walkers,1):
         xyz_grid_pts[rand_i[j]]=1
  if (Dim==2):
     for j in range(0,num_walkers,1):
         xyz_grid_pts[rand_i[j],rand_j[j]]=1
  if (Dim==3):
     for j in range(0,num_walkers,1):
         xyz_grid_pts[rand_i[j],rand_j[j],rand_k[j]]=1
  
  return xyz_grid_pts

def create_particles(Dim,xyz_grid_pts,num_walkers,rand_i,rand_j,rand_k,num_pts_x,num_pts_y,num_pts_z,probabilities):
  if (Dim==1):
     M=num_pts_x
     n=num_walkers
     for j in range(0,num_pts_x,1):
       if (np.random.uniform(0,1,1)<=(n/M)):
         xyz_grid_pts[j]=1
         n=n-1
         M=M-1
       else:
         xyz_grid_pts[j]=0
         M=M-1

  if (Dim==2):
     M=num_pts_x*num_pts_y
     n=num_walkers
     for i in range(0,num_pts_x,1):
       for j in range(0,num_pts_y,1):
         if (np.random.uniform(0,1,1)<=(n/M)):
           xyz_grid_pts[i,j]=1
           n=n-1
           M=M-1
         else:
           xyz_grid_pts[i,j]=0
           M=M-1


  if (Dim==3):
     M=num_pts_x*num_pts_y*num_pts_z
     n=num_walkers
     for i in range(0,num_pts_x,1):
       for j in range(0,num_pts_y,1):
         for k in range(0,num_pts_z,1):
           if (np.random.uniform(0,1,1)<=(n/M)):
             xyz_grid_pts[i,j,k]=1
             n=n-1
             M=M-1
           else:
             xyz_grid_pts[i,j,k]=0
             M=M-1

  return xyz_grid_pts


def move_particles(Dim,first_xyz_grid_pts,probabilities,timesteps,num_walkers,rand_i,rand_j,rand_k,num_pts_x,num_pts_y,num_pts_z,rand_indices):
  if (Dim==1):
    xyz_grid_pts=first_xyz_grid_pts
    count=[]
    for j in range(1,timesteps+1,1):
      particles_every_timestep=np.zeros(shape=(timesteps),dtype=int)
      rand_indices=np.asarray(np.where(xyz_grid_pts == 1))
      num_walkers=np.count_nonzero(xyz_grid_pts == 1)
      particle_num=np.zeros(shape=(num_walkers),dtype=int)
      particle_num=np.linspace(0,num_walkers-1,num_walkers)
      np.random.shuffle(particle_num)
      PN=particle_num.astype(int)
      num_walkers=np.count_nonzero(xyz_grid_pts == 1)
      for ip in range(0,num_walkers,1):
        rand_num=np.random.uniform(0,1,1)
        #print "random probability",rand_num
        #print "rand_indices\n",rand_indices
        #print "PN=",PN
        #print PN[ip], "numbered particle chosen"
        SY=rand_indices[0,PN[ip]]
        if (rand_num<=(1/2)):#move right
         #print "Move right probability = ",rand_num
         #print "(SY+1<num_pts_x) and (xyz_grid_pts[SY])",SY+1 , xyz_grid_pts[SY]
         if ((SY+1<num_pts_x) and (xyz_grid_pts[SY]!=0)):
          #print "xyz_grid_pts[SY+1]",xyz_grid_pts[SY+1]
          if ((xyz_grid_pts[SY+1]!=1)):
             #print "Entered move right loop and right slot is empty"
             xyz_grid_pts[SY]=0
             xyz_grid_pts[SY+1]=1
             #print "Right position updated"
          else:    
             #print "Entered move right loop and right slot is occupied"
             xyz_grid_pts[SY]=0
             xyz_grid_pts[SY+1]=0
             num_walkers=num_walkers-2
             #print "Indices of removed particles 1,2\n",SY,"\t",SY+1
             #print "Delete 2 particle right"
         else:
            pass
            #print "(SY+1<num_pts_x) and (xyz_grid_pts[SY]!=0) is FALSE"
            
        if (((1/2)<rand_num<=(2/2))):#move left
         #print "Move left probability = ",rand_num
         #print "(SY-1>=0) and (xyz_grid_pts[SY])",SY-1 , xyz_grid_pts[SY]
         if ((SY-1>=0) and (xyz_grid_pts[SY]!=0)): # and(xyz_grid_pts[SX,SY-1,SZ])
          #print "xyz_grid_pts[SY-1]",xyz_grid_pts[SY-1]
          if ((xyz_grid_pts[SY-1]!=1)):
             #print "Entered move left loop and left slot is empty"
             xyz_grid_pts[SY]=0
             xyz_grid_pts[SY-1]=1
             #print "Left position updated"
          else:
             #print "Entered move left loop and left slot is occupied"     
             xyz_grid_pts[SY]=0
             xyz_grid_pts[SY-1]=0
             #print "Delete 2 particles Left"
             num_walkers=num_walkers-2
             #print "Indices of removed particles 1,2\n",SY,"\t",SY-1
         else:
            pass
            #print "(SY-1>=0) and (xyz_grid_pts[SY]!=0) is FALSE"
            
      print np.count_nonzero(xyz_grid_pts == 1)  
      count.append(np.count_nonzero(xyz_grid_pts == 1))
      #particles_every_timestep[j-1]=np.count_nonzero(xyz_grid_pts == 1)
      #update_xyz[i-1,:]=xyz_grid_pts   
    particles_every_timestep=count

  if (Dim==2):
    count=[]
    xyz_grid_pts=first_xyz_grid_pts
    for j in range(1,timesteps+1,1):
      particles_every_timestep=np.zeros(shape=(timesteps,1),dtype=int)
      rand_indices=np.asarray(np.where(xyz_grid_pts == 1))
      num_walkers=np.count_nonzero(xyz_grid_pts == 1)
      particle_num=np.zeros(shape=(num_walkers),dtype=int)
      particle_num=np.linspace(0,num_walkers-1,num_walkers)
      np.random.shuffle(particle_num)
      PN=particle_num.astype(int)
      num_walkers=np.count_nonzero(xyz_grid_pts == 1)
      for ip in range(0,num_walkers,1):
        #print "After 1 particle moving, number of articles = ",np.count_nonzero(xyz_grid_pts == 1)
        rand_num=np.random.uniform(0,1,1)
        #print "random probability",rand_num
        #print "rand_indices = ",rand_indices
        #print "rand_indices[0,:],rand_indices[1,:]\n",rand_indices[0,:],rand_indices[1,:]
        #print "PN=",PN
        #print PN[ip], "numbered particle chosen"
        SX=rand_indices[0,PN[ip]]
        SY=rand_indices[1,PN[ip]]

        #print "SX,SY=",SX,SY
        if (rand_num<=(1/4)):#move right
         #print "Move right probability = ",rand_num
         #print "(SY+1<num_pts_y) and (xyz_grid_pts[SX,SY])",SY+1 , xyz_grid_pts[SX,SY]
         if ((SY+1<num_pts_y) and (xyz_grid_pts[SX,SY]!=0)):
          #print "xyz_grid_pts[SX,SY+1]",xyz_grid_pts[SX,SY+1]
          if ((xyz_grid_pts[SX,SY+1]!=1)):
             #print "Entered move right loop and right slot is empty"
             xyz_grid_pts[SX,SY]=0
             xyz_grid_pts[SX,SY+1]=1
             #print "Right position updated"
          else:    
             #print "Entered move right loop and right slot is occupied"
             xyz_grid_pts[SX,SY]=0
             xyz_grid_pts[SX,SY+1]=0
             num_walkers=num_walkers-2
             #print "Indices of removed particles 1,2\n",SX,SY,"\t",SX,SY+1
             #print "Delete 2 particle right"
         else:
            pass
            #print "(SY+1<num_pts_y) and (xyz_grid_pts[SX,SY]!=0) is FALSE"
            
        if (((1/4)<rand_num<=(2/4))):#move left
         #print "Move left probability = ",rand_num
         #print "(SY-1>=0) and (xyz_grid_pts[SX,SY])",SY-1 , xyz_grid_pts[SX,SY]
         if ((SY-1>=0) and (xyz_grid_pts[SX,SY]!=0)): 
          #print "xyz_grid_pts[SX,SY-1]",xyz_grid_pts[SX,SY-1]
          if ((xyz_grid_pts[SX,SY-1]!=1)):
             #print "Entered move left loop and left slot is empty"
             xyz_grid_pts[SX,SY]=0
             xyz_grid_pts[SX,SY-1]=1
             #print "Left position updated"
          else:
             #print "Entered move left loop and left slot is occupied"     
             xyz_grid_pts[SX,SY]=0
             xyz_grid_pts[SX,SY-1]=0
             #print "Delete 2 particles Left"
             num_walkers=num_walkers-2
             #print "Indices of removed particles 1,2\n",SX,SY,"\t",SX,SY-1
         else:
            pass
            #print "(SY-1>=0) and (xyz_grid_pts[SX,SY]!=0) is FALSE"
            
        if (((2/4)<rand_num<=(3/4))):#move down
         #print "Move down probability = ",rand_num
         #print "(SX+1<num_pts_x) and (xyz_grid_pts[SX,SY])",SX+1 , xyz_grid_pts[SX,SY]
         if ((SX+1<num_pts_x) and (xyz_grid_pts[SX,SY]!=0)): 
          #print "xyz_grid_pts[SX+1,SY]",xyz_grid_pts[SX+1,SY]
          if ((xyz_grid_pts[SX+1,SY]!=1)):
             #print "Entered move down loop and down slot is empty"
             xyz_grid_pts[SX,SY]=0
             xyz_grid_pts[SX+1,SY]=1
             #print "Down position updated"
          else:
             #print "Entered move down loop and down slot is occupied"        
             xyz_grid_pts[SX,SY]=0
             xyz_grid_pts[SX+1,SY]=0
             #print "Delete 2 particles down"
             num_walkers=num_walkers-2
             #print "Indices of removed particles 1,2\n",SX,SY,"\t",SX+1,SY
         else:
            pass
            #print "(SX+1<num_pts_x) and (xyz_grid_pts[SX,SY]!=0) is FALSE"
            

        if (((3/4)<rand_num<=(4/4))):#move up
         #print "Move up probability = ",rand_num
         #print "(SX-1>=0) and (xyz_grid_pts[SX,SY])",SX-1 , xyz_grid_pts[SX,SY]
         if ((SX-1>=0) and (xyz_grid_pts[SX,SY]!=0)):  
          #print "xyz_grid_pts[SX-1,SY]",xyz_grid_pts[SX-1,SY]
          if ((xyz_grid_pts[SX-1,SY]!=1)):
             #print "Entered move up loop and up slot is empty"
             xyz_grid_pts[SX,SY]=0
             xyz_grid_pts[SX-1,SY]=1
             #print "Up position updated"
          else:
             #print "Entered move up loop and up slot is occupied"          
             xyz_grid_pts[SX,SY]=0
             xyz_grid_pts[SX-1,SY]=0
             num_walkers=num_walkers-2
             #print "Indices of removed particles 1,2\n",SX,SY,"\t",SX-1,SY
             #print "Delete 2 particles down"
         else:
            pass
            #print "(SX-1<num_pts_x) and (xyz_grid_pts[SX,SY]!=0) is FALSE"
            

      print np.count_nonzero(xyz_grid_pts == 1)  
      count.append(np.count_nonzero(xyz_grid_pts == 1))
      #particles_every_timestep[j-1,0]=np.count_nonzero(xyz_grid_pts == 1)
      #update_xyz[i-1,:,:]=xyz_grid_pts      
    particles_every_timestep=np.asarray(count)

  if (Dim==3):
    count=[]
    xyz_grid_pts=first_xyz_grid_pts
    print "timestep 0", np.count_nonzero(xyz_grid_pts == 1)
    for j in range(1,timesteps+1,1):
      particles_every_timestep=np.zeros(shape=(timesteps),dtype=int)
      rand_indices=np.asarray(np.where(xyz_grid_pts == 1))
      #print "timestep",j
      num_walkers=np.count_nonzero(xyz_grid_pts == 1)
      particle_num=np.zeros(shape=(num_walkers),dtype=int)
      particle_num=np.linspace(0,num_walkers-1,num_walkers)
      np.random.shuffle(particle_num)
      PN=particle_num.astype(int)
      #print "PN=",PN
      num_walkers=np.count_nonzero(xyz_grid_pts == 1)
      for ip in range(0,num_walkers,1):
        #print "xyz_grid_pts=",xyz_grid_pts
        #print "After 1 particle moving, number of articles = ",np.count_nonzero(xyz_grid_pts == 1)
        rand_num=np.random.uniform(0,1,1)
        #print "random probability",rand_num
        #print "rand_indices = ",rand_indices
        #print "rand_indices[0,:],rand_indices[1,:],rand_indices[2,:]\n",rand_indices[0,:],rand_indices[1,:],rand_indices[2,:]
        #print "PN=",PN
        #print PN[ip], "numbered particle chosen"
        SX=rand_indices[0,PN[ip]]
        SY=rand_indices[1,PN[ip]]
        SZ=rand_indices[2,PN[ip]]

        #print "SX,SY,SZ=",SX,SY,SZ
        if (rand_num<=(1/6)):#move right
         #print "Move right probability = ",rand_num
         #print "(SY+1<num_pts_y) and (xyz_grid_pts[SX,SY,SZ])",SY+1 , xyz_grid_pts[SX,SY,SZ]
         if ((SY+1<num_pts_y) and (xyz_grid_pts[SX,SY,SZ]!=0)):
          #print "xyz_grid_pts[SX,SY+1,SZ]",xyz_grid_pts[SX,SY+1,SZ]
          if ((xyz_grid_pts[SX,SY+1,SZ]!=1)):
             #print "Entered move right loop and right slot is empty"
             xyz_grid_pts[SX,SY,SZ]=0
             xyz_grid_pts[SX,SY+1,SZ]=1
             #print "Right position updated"
             #print "Right move"
             #num_walkers=num_walkers
          else:    
#           if ((xyz_grid_pts[SX,SY+1,SZ]==1)):          
             #print "Entered move right loop and right slot is occupied"
             xyz_grid_pts[SX,SY,SZ]=0
             xyz_grid_pts[SX,SY+1,SZ]=0
             num_walkers=num_walkers-2
             #print "Indices of removed particles 1,2\n",SX,SY,SZ,"\t",SX,SY+1,SZ
             #print "Delete 2 particle right"
         else:
            pass
            #print "(SY+1<num_pts_y) and (xyz_grid_pts[SX,SY,SZ]!=0) is FALSE"
            
        if (((1/6)<rand_num<=(2/6))):#move left
         #print "Move left probability = ",rand_num
         #print "(SY-1>=0) and (xyz_grid_pts[SX,SY,SZ])",SY-1 , xyz_grid_pts[SX,SY,SZ]
         if ((SY-1>=0) and (xyz_grid_pts[SX,SY,SZ]!=0)): 
          #print "xyz_grid_pts[SX,SY-1,SZ]",xyz_grid_pts[SX,SY-1,SZ]
          if ((xyz_grid_pts[SX,SY-1,SZ]!=1)):
             #print "Entered move left loop and left slot is empty"
             xyz_grid_pts[SX,SY,SZ]=0
             xyz_grid_pts[SX,SY-1,SZ]=1
             #print "Left position updated"
             #print "Left move"
             #num_walkers=num_walkers
          else:   
             #print "Entered move left loop and left slot is occupied"     
             xyz_grid_pts[SX,SY,SZ]=0
             xyz_grid_pts[SX,SY-1,SZ]=0
             #print "Delete 2 particles Left"
             num_walkers=num_walkers-2
             #print "Indices of removed particles 1,2\n",SX,SY,SZ,"\t",SX,SY-1,SZ
         else:
            pass
            #print "(SY-1>=0) and (xyz_grid_pts[SX,SY,SZ]!=0) is FALSE"
            
        if (((2/6)<rand_num<=(3/6))):#move down
         #print "Move down probability = ",rand_num
         #print "(SX+1<num_pts_x) and (xyz_grid_pts[SX,SY,SZ])",SX+1 , xyz_grid_pts[SX,SY,SZ]
         if ((SX+1<num_pts_x) and (xyz_grid_pts[SX,SY,SZ]!=0)): 
          #print "xyz_grid_pts[SX+1,SY,SZ]",xyz_grid_pts[SX+1,SY,SZ]
          if ((xyz_grid_pts[SX+1,SY,SZ]!=1)):
             #print "Entered move down loop and down slot is empty"
             xyz_grid_pts[SX,SY,SZ]=0
             xyz_grid_pts[SX+1,SY,SZ]=1
             #print "Down position updated"
             #print "Down move"
             #num_walkers=num_walkers
          #if ((xyz_grid_pts[select_particle+1,select_particle,select_particle]!=1) and ((select_particle+1)>num_pts_x)):
             #xyz_grid_pts[select_particle,select_particle,select_particle]=1
          else:
          #if ((xyz_grid_pts[SX+1,SY,SZ]==1)):  
             #print "Entered move down loop and down slot is occupied"        
             xyz_grid_pts[SX,SY,SZ]=0
             xyz_grid_pts[SX+1,SY,SZ]=0
             #print "Delete 2 particles down"
             num_walkers=num_walkers-2
             #print "Indices of removed particles 1,2\n",SX,SY,SZ,"\t",SX+1,SY,SZ
         else:
            pass
            #print "(SX+1<num_pts_x) and (xyz_grid_pts[SX,SY,SZ]!=0) is FALSE"
            

        if (((3/6)<rand_num<=(4/6))):#move up
         #print "Move up probability = ",rand_num
         #print "(SX-1>=0) and (xyz_grid_pts[SX,SY,SZ])",SX-1 , xyz_grid_pts[SX,SY,SZ]
         if ((SX-1>=0) and (xyz_grid_pts[SX,SY,SZ]!=0)):  
          #print "xyz_grid_pts[SX-1,SY,SZ]",xyz_grid_pts[SX-1,SY,SZ]
          if ((xyz_grid_pts[SX-1,SY,SZ]!=1)):
             #print "Entered move up loop and up slot is empty"
             xyz_grid_pts[SX,SY,SZ]=0
             xyz_grid_pts[SX-1,SY,SZ]=1
             #print "Up position updated"
          else:
          #if ((xyz_grid_pts[SX-1,SY,SZ]==1)):
             #print "Entered move up loop and up slot is occupied"          
             xyz_grid_pts[SX,SY,SZ]=0
             xyz_grid_pts[SX-1,SY,SZ]=0
             num_walkers=num_walkers-2
             #print "Indices of removed particles 1,2\n",SX,SY,SZ,"\t",SX-1,SY,SZ
             #print "Delete 2 particles down"
         else:
            pass
            #print "(SX-1<num_pts_x) and (xyz_grid_pts[SX,SY,SZ]!=0) is FALSE"
            

        if (((4/6)<rand_num<=(5/6))):#move in
         #print "Move in probability = ",rand_num
         #print "(SZ-1>=0) and (xyz_grid_pts[SX,SY,SZ])",SZ-1 , xyz_grid_pts[SX,SY,SZ]
         if ((SZ-1>=0) and (xyz_grid_pts[SX,SY,SZ]!=0)): 
          #print "xyz_grid_pts[SX,SY,SZ-1]",xyz_grid_pts[SX,SY,SZ-1]
          if ((xyz_grid_pts[SX,SY,SZ-1]!=1)):
             #print "Entered move in loop and in slot is empty"
             xyz_grid_pts[SX,SY,SZ]=0
             xyz_grid_pts[SX,SY,SZ-1]=1
             #print "In position updated"
             #print "In move"
             #num_walkers=num_walkers

          else:
          #if ((xyz_grid_pts[SX,SY,SZ-1]==1)):  
             #print "Entered move in loop and in slot is occupied"        
             xyz_grid_pts[SX,SY,SZ]=0
             xyz_grid_pts[SX,SY,SZ-1]=0
             #print "Delete 2 particles in"
             num_walkers=num_walkers-2
             #print "Indices of removed particles 1,2\n",SX,SY,SZ,"\t",SX,SY,SZ-1
         else:
            pass
            #print "(SZ-1<num_pts_z) and (xyz_grid_pts[SX,SY,SZ]!=0) is FALSE"
            

        if (((5/6)<rand_num<=(6/6))):#move out
         #print "Move out probability = ",rand_num
         #print "(SZ+1<num_pts_z) and (xyz_grid_pts[SX,SY,SZ])",SZ+1 , xyz_grid_pts[SX,SY,SZ]
         if ((SZ+1<num_pts_z) and (xyz_grid_pts[SX,SY,SZ]!=0)):
          #print "xyz_grid_pts[SX,SY,SZ+1]",xyz_grid_pts[SX,SY,SZ+1]
          if ((xyz_grid_pts[SX,SY,SZ+1]!=1)):
             #print "Entered move out loop and out slot is empty"
             xyz_grid_pts[SX,SY,SZ]=0
             xyz_grid_pts[SX,SY,SZ+1]=1
             #print "Out position updated"
          else:
          #if ((xyz_grid_pts[SX,SY,SZ+1]==1)):   
             #print "Entered move out loop and out slot is occupied"       
             xyz_grid_pts[SX,SY,SZ]=0
             xyz_grid_pts[SX,SY,SZ+1]=0
             num_walkers=num_walkers-2
             #print "Delete 2 particles out"
             #print "Indices of removed particles 1,2\n",SX,SY,SZ,"\t",SX,SY,SZ-1
         else:
            pass
            #print "(SZ+1<num_pts_z) and (xyz_grid_pts[SX,SY,SZ]!=0) is FALSE"
            

        
        #print "xyz_grid_pts\n",xyz_grid_pts
      print np.count_nonzero(xyz_grid_pts == 1)  
      count.append(np.count_nonzero(xyz_grid_pts == 1))
    particles_every_timestep=np.asarray(count)    

  return xyz_grid_pts,particles_every_timestep#update_xyz

timesteps=100000

#xyz_grid_pts=create_random_particles(Dim,xyz_grid_pts,num_walkers,rand_i,rand_j,rand_k,probabilities)
first_xyz_grid_pts=create_particles(Dim,xyz_grid_pts,num_walkers,rand_i,rand_j,rand_k,num_pts_x,num_pts_y,num_pts_z,probabilities)


num_ensembles=600

if (Dim==1):
 M=(num_pts_x)  ### Number of lattice sites

if (Dim==2):
 M=(num_pts_x*num_pts_y)  ### Number of lattice sites

if (Dim==3):
 M=(num_pts_x*num_pts_y*num_pts_z)  ### Number of lattice sites

#print num_ensembles
#print first_xyz_grid_pts  

particles=np.count_nonzero(xyz_grid_pts == 1)
#print particles
rand_indices=np.asarray(np.where(xyz_grid_pts == 1))
size=np.asarray(rand_indices.shape)

update_xyz,particles_every_timestep=move_particles(Dim,first_xyz_grid_pts,probabilities,timesteps,num_walkers,rand_i,rand_j,rand_k,num_pts_x,num_pts_y,num_pts_z,rand_indices)

np.savetxt('Number_of_particles_at_every_timestep.dat',particles_every_timestep)

Ens_atoms_at_t=np.zeros(shape=(timesteps),dtype=float)
for E in range(1,num_ensembles+1,1):   
   print "Ensemble", E
   first_xyz_grid_pts=create_particles(Dim,xyz_grid_pts,num_walkers,rand_i,rand_j,rand_k,num_pts_x,num_pts_y,num_pts_z,probabilities)

   update_xyz,particles_every_timestep=move_particles(Dim,first_xyz_grid_pts,probabilities,timesteps,num_walkers,rand_i,rand_j,rand_k,num_pts_x,num_pts_y,num_pts_z,rand_indices)

   Ens_atoms_at_t=Ens_atoms_at_t+particles_every_timestep

Ens_dens_avg=np.zeros(shape=(timesteps),dtype=float)
Ens_dens_avg=Ens_atoms_at_t/(M*num_ensembles)   

np.savetxt('Ensemble_density_average_at_every_timestep.dat',Ens_dens_avg)

if (Dim==1):
  data=np.loadtxt('Ensemble_density_average_at_every_timestep.dat') #Ensemble average density values at different times
  density=data
  timestep=np.linspace(1,timesteps,timesteps) ## The values depend on chosen number of timesteps

###To plot analytic curves
  x=np.linspace(1,1001,1001)
  y=np.zeros(shape=(1001),dtype=float)
  for i in range(0,1001,1):
   y[i]=1/sqrt(x[i])

  plot(Ens_dens_avg,timestep,'r-o',label='Simulation, 1D Lattice size = 100000 units\n800 walkers at t=0')
  #plot(x,y,'y-x',label='Analytical')
  ylabel('Density of walkers') 
  xlabel('Timestep')
  legend(loc='best')

  savefig('1D_particles_timestep.png')
  #savefig('1D_Analytical_curve.png')

if (Dim==2):
  data=np.loadtxt('Ensemble_density_average_at_every_timestep.dat') #Ensemble average density values at different times
  density=data
  timestep=np.linspace(1,timesteps,timesteps) ## The values depend on chosen number of timesteps

###To plot analytic curves
  x=np.linspace(1,1001,1001)
  y=np.zeros(shape=(1001),dtype=float)
  for i in range(0,1001,1):
    y[i]=log(x[i])/x[i]

  plot(Ens_dens_avg,timestep,'r-o',label='Simulation, 2D Lattice size = 500X500\n5000 walkers at t=0')
  #plot(x,y,'y-x',label='Analytical')
  ylabel('Density of walkers') 
  xlabel('Timestep')
  legend(loc='best')

  savefig('2D_particles_timestep.png')
  #savefig('2D_Analytical_curve.png')

if (Dim==3):
  data=np.loadtxt('Ensemble_density_average_at_every_timestep.dat') #Ensemble average density values at different times
  density=data
  timestep=np.linspace(1,timesteps,timesteps) ## The values depend on chosen number of timesteps

###To plot analytic curves
  x=np.linspace(1,1001,1001)
  y=np.zeros(shape=(1001),dtype=float)
  for i in range(0,1001,1):
     y[i]=1/x[i]

  plot(Ens_dens_avg,timestep,'r-o',label='Simulation, 3D Lattice size = 200X200X200\n343 walkers at t=0')
  #plot(x,y,'y-x',label='Analytical')
  ylabel('Density of walkers') 
  xlabel('Timestep')
  legend(loc='best')

  savefig('3D_particles_timestep.png')
  #savefig('3D_Analytical_curve.png')





