import numpy as np 

cdef Acceleration(coordFirst, coordSecond, mass):
	cdef double G = 6.67408e-11
	if (np.linalg.norm(coordSecond - coordFirst) < 0.0000001):
		return 0
	return G*mass*(coordSecond - coordFirst)/(np.linalg.norm(coordSecond - coordFirst)**3)

cdef VerletStep(R, V, A, i, dt, kwargs):
	for (mass, r) in kwargs:
		A[:, i] += Acceleration(R[:, i], r[:, i], mass)
	R[:, i+1] = R[:, i] + V[:, i]*dt + (A[:, i]*(dt**2)*0.5)
	V[:, i+1] = V[:, i] + A[:, i]*dt

class Verlet:
	def __init__(self, Num, Dt, objects):
		self.N = Num
		self.dt = Dt
		self.Objects = objects
	def VerletMain(self):
		for i in range(self.N):
			for obj in self.Objects:
				VerletStep(obj.R, obj.V, obj.A, i, self.dt, obj.Interactions)

cdef VerletStepOpenMP(R, V, A, i, dt, kwargs):
	#pragma parallel for
	for (mass, r) in kwargs:
		A[:, i] += Acceleration(R[:, i], r[:, i], mass)
	R[:, i+1] = R[:, i] + V[:, i]*dt + (A[:, i]*(dt**2)*0.5)
	V[:, i+1] = V[:, i] + A[:, i]*dt

class VerletOpenMP:
	def __init__(self, Num, Dt, objects):
		self.N = Num
		self.dt = Dt
		self.Objects = objects
	def VerletMain(self):
		for i in range(self.N):
			#pragma parallel for
			for obj in self.Objects:
				VerletStepOpenMP(obj.R, obj.V, obj.A, i, self.dt, obj.Interactions)

cdef double Norm(double[:] vector1, double[:] vector2):
	cdef double res = 0
	cdef int i
	for i in range(3):
		res += (vector1[i]-vector2[i])**2
	res = res ** 0.5
	return res

cdef double[:] AccelerationTM(double[:] coordFirst, double[:] coordSecond,double mass):
	cdef double G = 6.67408e-11
	cdef int i
	cdef double[:] res = np.zeros(3)
	if (Norm(coordFirst, coordSecond) < 0.0000001):
		return res
	for i in range(3):
		res[i] = G*mass*(coordSecond[i] - coordFirst[i])/(Norm(coordFirst, coordSecond)**3)
	return res

cdef VerletStepTM(double[:,:] R, double[:,:] V, double[:,:] A, int i, int dt, kwargs):
	cdef int k
	cdef double[:] APart
	for (mass, r) in kwargs:
		APart = AccelerationTM(R[:, i], r[:, i], mass)
		for k in range (3):
			A[k, i] += APart[k]
	for k in range (3):
		R[k, i+1] = R[k, i] + V[k, i]*dt + (A[k, i]*(dt**2)*0.5)
		V[k, i+1] = V[k, i] + A[k, i]*dt

class VerletTM:
	def __init__(self, Num, Dt, objects):
		self.N = Num
		self.dt = Dt
		self.Objects = objects
	def VerletMain(self):
		for i in range(self.N):
			for obj in self.Objects:
				VerletStepTM(obj.R, obj.V, obj.A, i, self.dt, obj.Interactions)

cdef VerletStepOpenMPTM(double[:,:] R, double[:,:] V, double[:,:] A, int i, int dt, kwargs):
	cdef int k
	cdef double[:] APart
	#pragma parallel for
	for (mass, r) in kwargs:
		APart = AccelerationTM(R[:, i], r[:, i], mass)
		for k in range (3):
			A[k, i] += APart[k]
	for k in range (3):
		R[k, i+1] = R[k, i] + V[k, i]*dt + (A[k, i]*(dt**2)*0.5)
		V[k, i+1] = V[k, i] + A[k, i]*dt

class VerletOpenMPTM:
	def __init__(self, Num, Dt, objects):
		self.N = Num
		self.dt = Dt
		self.Objects = objects
	def VerletMain(self):
		for i in range(self.N):
			#pragma parallel for
			for obj in self.Objects:
				VerletStepOpenMPTM(obj.R, obj.V, obj.A, i, self.dt, obj.Interactions)