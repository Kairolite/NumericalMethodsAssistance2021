import numpy as np
from matplotlib import pyplot as plt

#Data variance
X = np.array([0,2,4,6,9,11,12,15,17,19])
Y = np.array([5,6,7,6,9,8,8,10,12,12])

#Regression order
Order = 3
n = Order+1

#Matrix of z
Z = np.zeros(shape=(n, X.size))
for i in reversed(range(n)):
    Z[i]= X**(Order-i)
Ztr = np.transpose(Z)
A = Z.dot(Ztr)
B = Z.dot(Y) 

#Creating Augmented Matrix
a = np.zeros((n,n+1))

for i in range(n):
    a[:,i] = A[:,i]
a[:,-1]=B

#print (Z)
#print (Ztr)
#print ("A =")
print (A)
#print ("B =")
print ("")
print (B)
print ("")

print ("Gauss Jordan =")
print (a)
print ("")

#Gauss Jordan Elimination, only trivial matrix
x = np.zeros(n)
for i in range(n):
    for j in range(n):
        if i != j:
            ratio = a[j][i]/a[i][i]

            for k in range(n+1):
            	a[j][k] = a[j][k] - ratio * a[i][k]
            	
for i in range(n):
	x[i] = a[i][n]/a[i][i]

print('The solutions are: ')

for i in range(n):
    print('X%d = %0.5f' %(i+1,x[i]), end = ' ')

x_ = np.linspace(0, 20, 21)
y_ = x[0]*x_**3 + x[1]*x_**2 + x[2]*x_**1 + x[3]


plt.scatter(X, Y)
#plt.plot(X, Y)
plt.plot(x_, y_, '--', color="red")
plt.show()
