#include "colors.inc"
#include "woods.inc"
#include "stones.inc"
#include "metals.inc"
#include "golds.inc"
#include "glass.inc"

//global_settings{ assumed_gamma 1.0 }
//#default{ finish{ ambient 0.1 diffuse 0.9 }} 

//Place the camera
camera {
  sky <0,0,1>          
  location  <0, 2, -3>  
  look_at   <0.0,0.4,0.0>   
}

//Place a light
light_source {
  <0, 10,-10>  //Change this if you want to put the light at a different point
  color White
}

light_source {
  <0, 0,-100>  //Change this if you want to put the light at a different point
  color 0.5*White
}

#declare v0 = < -1,0,-1>;
#declare v1 = <  1,0,-1>;
#declare v2 = <  1,0, 1>;
#declare v3 = < -1,0, 1>;
#declare v4 = < -0,1.5, 0>;
#declare v5 = 0.5*(v0 + v1 );
#declare v6 = 0.5*(v1 + v2 );
#declare v7 = 0.5*(v2 + v3 );
#declare v8 = 0.5*(v3 + v0 );

//Create a "floor"
plane {
  <0,1,0>,-1.1
  texture { T_Stone5 }
}

#declare cylRadius = 0.05;
#declare sphRadius = 0.1;

#declare node0 = sphere{v0, sphRadius};
#declare node1 = sphere{v1, sphRadius};
#declare node2 = sphere{v2, sphRadius};
#declare node3 = sphere{v3, sphRadius};
#declare node4 = sphere{v4, sphRadius};
#declare node5 = sphere{v5, sphRadius};
#declare node6 = sphere{v6, sphRadius};
#declare node7 = sphere{v7, sphRadius};
#declare node8 = sphere{v8, sphRadius};

#declare edge1 = cylinder{v0,v1, cylRadius};
#declare edge2 = cylinder{v1,v2, cylRadius};
#declare edge3 = cylinder{v2,v3, cylRadius};
#declare edge4 = cylinder{v3,v0, cylRadius};

#declare edge5 = cylinder{v0,v4, cylRadius};
#declare edge6 = cylinder{v1,v4, cylRadius};
#declare edge7 = cylinder{v2,v4, cylRadius};
#declare edge8 = cylinder{v3,v4, cylRadius};

#declare HexNodes = 
union { 
       object{node0}
       object{node1}
       object{node2}
       object{node3}
       object{node4}
       object{node5}
       object{node6}
       object{node7}
       object{node8}
       texture {T_Gold_1A}
}

#declare HexEdges = 
union { 
       object{edge1}
       object{edge2}
       object{edge3}
       object{edge4}
       object{edge5}
       object{edge6}
       object{edge7}
       object{edge8}
       texture {T_Silver_1A}
}
  
#declare Scene = 
union {
     object {HexNodes}
     object {HexEdges}
     rotate<0, 360*clock, 0>
}

Scene
