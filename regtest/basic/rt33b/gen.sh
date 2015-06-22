awk 'BEGIN{
for(i=0;i<100;i++){
  print "3"
  print "0 0 0"
  print "Ar 0 0 0"
  print "Ar ",rand(), rand(),rand();
  print "Ar ",rand(), rand(),rand();
}
}' > trajectory.xyz

{
cat << EOF
#! FIELDS time d.x sigma_d.x height biasf
#! SET multivariate false
EOF
awk 'BEGIN{
for(i=0;i<100;i++){
  print i,rand(),(0.05+rand()*0.1),1.0,1.0
}
}'
} > HILLS1

awk '{
  if(match($1,"#"))print
  else print $1,$2,$3,-$4,$5
}' HILLS1 > HILLS1m

{
cat << EOF
#! FIELDS time d.y d.z sigma_d.y sigma_d.z height biasf
#! SET multivariate false
EOF
awk 'BEGIN{
for(i=0;i<100;i++){
  print i,rand(),rand(),(0.05+rand()*0.1),(0.05+rand()*0.1),1.0,1.0
}
}'
} > HILLS2

awk '{
  if(match($1,"#"))print
  else print $1,$2,$3,$4,$5,-$6,$7
}' HILLS2 > HILLS2m

{
cat << EOF
#! FIELDS time c.x c.y c.z sigma_c.x sigma_c.y sigma_c.z height biasf
#! SET multivariate false
EOF
awk 'BEGIN{
for(i=0;i<100;i++){
  print i,rand(),rand(),rand(),(0.05+rand()*0.1),(0.05+rand()*0.1),(0.05+rand()*0.1),1.0,1.0
}
}'
} > HILLS3

awk '{
  if(match($1,"#"))print
  else print $1,$2,$3,$4,$5,$6,$7,-$8,$9
}' HILLS3 > HILLS3m



