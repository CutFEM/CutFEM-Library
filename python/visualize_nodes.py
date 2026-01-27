# import numpy as np
# import matplotlib.pyplot as plt

# # Files produced by your code
# vol_file  = "../build/cuttri_multipoly_vol_quad.txt"
# surf_file = "../build/cuttri_multipoly_surf_quad.txt"

# vol  = np.loadtxt(vol_file)   # columns: x y w
# surf = np.loadtxt(surf_file)

# xv, yv, wv = vol[:,0],  vol[:,1],  vol[:,2]
# xs, ys, ws = surf[:,0], surf[:,1], surf[:,2]

# plt.figure()
# plt.gca().set_aspect("equal", adjustable="box")

# # Volume points
# # plt.scatter(xv, yv, s=30*(wv/np.max(wv)), alpha=0.6, label="volume nodes")
# plt.scatter(xv, yv, s=10, alpha=0.6, label="volume nodes")

# # Surface points
# plt.scatter(xs, ys, s=60*(ws/np.max(ws)), alpha=0.9, label="surface nodes")

# plt.legend()
# plt.title("Algoim multipoly quadrature nodes (size ∝ weight)")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.grid(True)
# plt.show()



import numpy as np
import matplotlib.pyplot as plt

vol  = np.loadtxt("../build/cuttri_multipoly_vol_quad.txt")
surf = np.loadtxt("../build/cuttri_multipoly_surf_quad.txt")

xv, yv, wv = vol[:,0],  vol[:,1],  vol[:,2]
xs, ys, ws = surf[:,0], surf[:,1], surf[:,2]

plt.figure()
plt.gca().set_aspect("equal", adjustable="box")

# triangle
# tri = np.array([[0,0],[1,0],[0,1],[0,0]])
tri = np.array([[0,0],[0.33,0.33],[0,1],[0,0]])
plt.plot(tri[:,0], tri[:,1], "-k", lw=2)

# circle boundary (center (0.5,0.5), r=0.6)
t = np.linspace(0, 2*np.pi, 400)
plt.plot(0.5 + 0.6*np.cos(t), 0.5 + 0.6*np.sin(t), "--k", lw=1)

plt.scatter(xv, yv, s=8, label="volume nodes")
plt.scatter(xs, ys, s=25, label="surface nodes")
plt.legend()
plt.grid(True)
plt.show()

print("max(x+y) =", (xv+yv).max())     # should be <= 1 (up to tiny tol) for triangle (0,0)-(1,0)-(0,1)
print("sum(w)   =", wv.sum())         # area approximation of K ∩ {phi<0}
print("min(w)   =", wv.min(), "max(w) =", wv.max())

