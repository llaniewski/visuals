#setwd("~/ludziki")
rm(list=ls())

bigdim_nx = 300
bigdim_ny = bigdim_nx
tiles_nx = 1
tiles_ny = 1
tile = 1
subsamp = 5
subsamp_p = 4 # power
N = 7
M = 3
tp = "idn"


move = function(z,w) (z-w)/(1-w*z)
rota = function(z,al) z*exp(1i*al)
flip = function(z) -Re(z)+1i*Im(z)

nicediv = function(a,b) {
  m = a %% b
  n = rep((a-m)/b,b)
  n[seq_len(m)] = n[seq_len(m)] + 1
  n
}

alpha = 2*pi/N
beta = 2*pi/M
A = cos(alpha/2)/sin(beta/2)
W = sqrt((A-1)/(A+1))

p = expand.grid(tile_y=1:tiles_ny,tile_x=1:tiles_nx)
a = nicediv(bigdim_nx, tiles_nx)
p$tile_nx = a[p$tile_x]
p$tile_dx = c(0,cumsum(a))[p$tile_x]
a = nicediv(bigdim_ny, tiles_ny)
p$tile_ny = a[p$tile_y]
p$tile_dy = c(0,cumsum(a))[p$tile_y]

for (tile in seq_len(nrow(p))) {

if (tile > nrow(p)) stop("tile number out of range")

dim_nx = p$tile_nx[tile]
dim_ny = p$tile_ny[tile]
dim_dx = p$tile_dx[tile]
dim_dy = p$tile_dy[tile]

library(png)

hips = function(z, an) {
  z = rota(z,an)
  z = move(z,W)
  z = flip(z)
  z = move(z,-W)
  z = rota(z,-an)
  z
}

x = seq_len(bigdim_nx)  - (bigdim_nx+1)/2
y = seq_len(bigdim_ny)  - (bigdim_ny+1)/2
scale = 2 / bigdim_ny
x = x * scale
y = y * scale
x = x[seq_len(dim_nx)+dim_dx]
y = y[seq_len(dim_ny)+dim_dy]
z = expand.grid(x=x,y=y)
z = z$x + z$y * 1i

k = subsamp
dz = expand.grid(1:k,1:k)-0.5
dz = (dz[,1]+1i*dz[,2])/k
#dx = runif(k*k) + runif(k*k)*1i
dz = runif(k*k) + runif(k*k)*1i
dz = dz - (0.5 + 0.5i)
dz = dz * scale
subsamp_n = length(dz)

z = expand.grid(z=z,dz=dz)
z = z$z + z$dz

z = z*1.1

if (tp == "hip") {
f = z

pminm = function(x,y) ifelse(Mod(x) < Mod(y), x, y)
minm = function(x) {
  r = x[,1]
  for (i in 2:ncol(x)) r = pminm(x[,i],r)
  r
}

self=seq_along(f)

self = self[Mod(f)<1]

for (i in 1:30) {
  print(length(self))
  sf = f[self]
  nf = cbind(sf,do.call(cbind, lapply(1:M, function(i) hips(sf,2*pi/M*i))))
  nf = minm(nf)
  f[self] = nf
  self = self[nf != sf]
  if (length(self) <= 0) break;
}

getw = function(f) {
  x = Re(f)
  y = Im(f)
  A = 1+x^2+y^2
  2*x/(A+sqrt(A*A-4*x*x))
}

self = FALSE
for (i in 1:M) self = self | getw(rota(f,2*pi/M*i)) > W*0.9
f = getw(f) + getw(rota(f,pi/2))*1i
f = f/W/2
lum = rep(0.9,dim_nx*dim_ny*subsamp_n)
lum[self] = 0

} else if (tp == "inv") {
f = rota(1/z,0.2)
} else if (tp == "idn") {
  p = cbind(Re(z),Im(z))
  anx = 0.3
  np = p*9.888
  np = np %*% matrix(c(cos(anx),sin(anx),-sin(anx),cos(anx)),2,2)
  W = cbind(cos(2*pi*(1:3)/3),sin(2*pi*(1:3)/3))
  rad = function(p) p[,1]^2+p[,2]^2+(p[,3]-1)^2
  for (j in 1:20) {
    for (i in 1:nrow(W)) {
      v = W[i,]
      vp = np - 2*(np %*% v-1) %*% v
      sel = rowSums(np^2) > rowSums(vp^2)
      print(mean(sel))
      np[sel] = vp[sel]
    }
  }
  f = np[,1]+1i*np[,2]
  sel = FALSE
  for (i in 1:nrow(W)) {
    v = W[i,]
    sel = sel | (1-np %*% v) < 0.1
  }
  lum = rep(0.9,dim_nx*dim_ny*subsamp_n)
  lum[sel] = 0
} else if (tp == "ico") {
  q = z[Mod(z) < 1]
  p = cbind(Re(q),Im(q), sqrt(1-Mod(q)^2))
  anx = 0.3
  any = 0.3
  np = p
  np = np %*% matrix(c(cos(any),0,sin(any),0,1,0,-sin(any),0,cos(any)),3,3)
  np = np %*% matrix(c(1,0,0,0,cos(anx),sin(anx),0,-sin(anx),cos(anx)),3,3)
  W = { a = 2*pi/5
        b = (1:3-1)*2*pi/3
        s = 1/sin(2*pi/5)*sin(2*pi/6)*2/3
        c = sqrt(1-s^2)
        V = cbind(s*cos(b),s*sin(b),c)
        W = t(solve(V))
        W / sqrt(rowSums(W^2))
  }
  rad = function(p) p[,1]^2+p[,2]^2+(p[,3]-1)^2
  for (j in 1:20) {
    for (i in 1:nrow(W)) {
      v = W[i,]
      vp = np - 2*(np %*% v) %*% v
      sel = rad(np) > rad(vp)
      print(mean(sel))
      np[sel] = vp[sel]
    }
  }
  f = z
  f[Mod(z) < 1] = np[,1]+1i*np[,2]
  sel = FALSE
  for (i in 1:nrow(W)) {
    v = W[i,]
    sel = sel | np %*% v < 0.03
  }
  L = c(-1,1,1)
  L = L/sqrt(sum(L^2))
  L = outer(rep(1,nrow(p)),L)
  N = p
  R = 2 * rowSums(N * L) * N - L
  V = outer(rep(1,nrow(p)),c(0,0,1))
  lum = rep(1,dim_nx*dim_ny*subsamp_n)
  i1 = rep(0.6,nrow(p))
  i2 = rep(0.3,nrow(p))
  i3 = rep(0.1,nrow(p))
  i1[sel] = 0.2
  i2[sel] = 0
  lum[Mod(z) < 1] = i1 + i2*pmax(0,rowSums(N * L)) + i3*pmax(0,rowSums(R * V))^4
#  lum[Mod(z) < 1][sel] = 0
} else {
  stop("unknown")
}

x = Re(f)-0.5
y = Im(f)-0.5

symx = floor(x) %% 2 == 0
symy = floor(y) %% 2 == 0

x = x %% 1
y = y %% 1

x = ifelse(symx,x,1-x)
y = ifelse(symy,y,1-y)

input = readPNG("c_test.png")
if (length(dim(input)) < 3) dim(input) = c(dim(input),1)
input_nx = dim(input)[1]
input_ny = dim(input)[2]
channels = dim(input)[3]
dim(input) = c(input_nx*input_nx, channels)
output = matrix(0, dim_nx*dim_ny, channels)

for (channel in seq_len(channels)) {
#  h = input[floor(x * input_nx) + floor(y * input_ny) * input_nx + 1, channel]
  h = lum
  h[Mod(z) > 1] = 1
  dim(h) = c(dim_nx*dim_ny, subsamp_n)
  h = 1-rowMeans((1-h)^subsamp_p)^(1/subsamp_p)
  output[,channel] = h
}
if (channels == 1) {
  dim(output) = c(dim_nx,dim_ny)
} else {
  dim(output) = c(dim_nx,dim_ny,channels)
}

writePNG(output,target=sprintf("img_%s_%03d.png",tp,tile))
}





