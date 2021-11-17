import numpy as np
from scipy import io, integrate, linalg, signal
from scipy.sparse.linalg import eigs
import matplotlib.pyplot as plt
import matplotlib.image as img
#camera matrix
K=np.array([[-3933.36363636,0,2143.5],[0,3933.36363636,1423.5],[0,0,1]])
#roation matrix
R=np.array([[-0.827570749773,-0.557530411593,0.065471324026],[0.549857636160,-0.828569770064,-0.105492730048],[0.113062965097,-0.051302790236,0.992262460056]])
#translation vector
t=np.array([[512980.995110],[5427701.526710],[514.794290]])

#projection matrix
P=np.dot(K,np.c_[R,-np.dot(R,t)])

#object coordinate
X=np.array([[512868.940,5427723.833,280.885],[512982.899,5427801.281,323.376],[513047.601,5427704.443,326.998],[512781.152,5427721.195,229.286],[512866.691,5427556.801,229.433],[512990.398,5427681.815,330.030],[512995.967,5427683.424,326.758]])
#homogen coordiante
X=np.c_[X,np.ones((7,1))]

x=np.dot(P,X.T)
#nomination
x_n=np.vstack((x[0,:]/x[2,:],x[1,:]/x[2,:],np.ones(7).T)).T

#load image
# image=img.imread('R0020851.jpg')
# plt.imshow(image)
# plt.scatter(x_n[:,0],x_n[:,1])
# plt.show()

#task 2
A=np.zeros((14,12))
for i in range(1,8):
	A[2*i-2,:]=np.hstack((np.zeros(4),-X[i-1,:],x_n[i-1,1]*X[i-1,:]))
	A[2*i-1,:]=np.hstack((X[i-1,:],np.zeros(4),-x_n[i-1,0]*X[i-1,:]))
# print(A)
U, S, Vh = linalg.svd(A,0) 
V = Vh.T
P2= V[:,11].reshape(4,3,order='F').copy().T

x2=np.dot(P2,X.T)
x2=np.vstack((x2[0,:]/x2[2,:],x2[1,:]/x2[2,:])).T
# image=img.imread('R0020851.jpg')
# plt.imshow(image)
# plt.scatter(x2[:,0],x2[:,1])
# plt.show()

#task 3
M=P[0:3,0:3]
[RotTrans,Kinv]=linalg.qr(linalg.inv(M))
R_re=RotTrans.T
K_re=linalg.inv(Kinv)
K_re=K_re/K_re[2,2]

print(R_re)
print(K_re)
