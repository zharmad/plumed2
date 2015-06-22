awk '{
  if(match($0,"#"))print;
  else {
    x=$1
    print x,x**3-x**2+2,3*x**2-2*x
  }
}' grid1 | sed 's/e1/e1p/' > grid1poly

awk '{
  if(match($0,"#"))print;
  else if (NF==0) print 
  else {
    x=$1
    y=$2
    print x,y, x**3-x**2 + y**2, 3*x**2-2*x, 2*y
  }
}' grid2 | sed 's/e2/e2p/' > grid2poly

awk '{
  if(match($0,"#"))print;
  else if (NF==0) print 
  else {
    x=$1
    y=$2
    z=$3
    print x,y,z,  x+y**3+y-z**3, 1 , 3*y**2 +1 , -3*z**2
  }
}' grid3 | sed 's/e3/e3p/' > grid3poly


