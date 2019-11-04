#include "colors.inc"
#include "woods.inc"
#include "stones.inc"
#include "metals.inc"
#include "golds.inc"
#include "glass.inc"

camera {
  sky <0,0,1>          
  location  <0,  1, -5>  
  look_at   <0.0,1.0,0.0>   
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

#declare len = 1.0;

#declare v0 = <-len,0,len>;
#declare v1 = < len,0,len>;
#declare v2 = < len,2*len,len>;
#declare v3 = < -len,2*len,len>;

#declare v4 = <-len,0,-len>;
#declare v5 = < len,0,-len>;
#declare v6 = < len,2*len,-len>;
#declare v7 = < -len,2*len,-len>;

#declare v8  = <- 0.5*len, 0.5*len, 0.5*len>;
#declare v9  = <  0.5*len, 0.5*len, 0.5*len>;
#declare v10 = <  0.5*len, 1.5*len, 0.5*len>;
#declare v11 = < -0.5*len, 1.5*len, 0.5*len>;

#declare v12 = <- 0.5*len, 0.5*len, -0.5*len>;
#declare v13 = <  0.5*len, 0.5*len, -0.5*len>;
#declare v14 = <  0.5*len, 1.5*len, -0.5*len>;
#declare v15 = < -0.5*len, 1.5*len, -0.5*len>;

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

#declare node8  = sphere{v8, sphRadius};
#declare node9  = sphere{v9, sphRadius};
#declare node10 = sphere{v10, sphRadius};
#declare node11 = sphere{v11, sphRadius};
#declare node12 = sphere{v12, sphRadius};
#declare node13 = sphere{v13, sphRadius};
#declare node14 = sphere{v14, sphRadius};
#declare node15 = sphere{v15, sphRadius};


// Big Cube
#declare edge1 = cylinder{v0,v1, cylRadius};
#declare edge2 = cylinder{v0,v3, cylRadius};
#declare edge3 = cylinder{v0,v4, cylRadius};
#declare edge4 = cylinder{v1,v2, cylRadius};
#declare edge5 = cylinder{v1,v5, cylRadius};
#declare edge6 = cylinder{v2,v3, cylRadius};
#declare edge7 = cylinder{v2,v6, cylRadius};
#declare edge8 = cylinder{v3,v7, cylRadius};
#declare edge9  = cylinder{v4,v5, cylRadius};
#declare edge10 = cylinder{v4,v7, cylRadius};
#declare edge11 = cylinder{v5,v6, cylRadius};
#declare edge12 = cylinder{v6,v7, cylRadius};

// Bridges ...
#declare edge13 = cylinder{v0,v8, cylRadius};
#declare edge14 = cylinder{v1,v9, cylRadius};
#declare edge15 = cylinder{v2,v10, cylRadius};
#declare edge16 = cylinder{v3,v11, cylRadius};

#declare edge17 = cylinder{v4,v12, cylRadius};
#declare edge18 = cylinder{v5,v13, cylRadius};
#declare edge19 = cylinder{v6,v14, cylRadius};
#declare edge20 = cylinder{v7,v15, cylRadius};

// Smaller Cube

#declare edge21 = cylinder{v8,v9, cylRadius};
#declare edge22 = cylinder{v8,v11, cylRadius};
#declare edge23 = cylinder{v8,v12, cylRadius};

#declare edge24 = cylinder{v9,v10, cylRadius};
#declare edge25 = cylinder{v9,v13, cylRadius};

#declare edge26 = cylinder{v10,v11, cylRadius};
#declare edge27 = cylinder{v10,v14, cylRadius};

#declare edge28 = cylinder{v11,v15, cylRadius};

#declare edge29 = cylinder{v12,v13, cylRadius};
#declare edge30 = cylinder{v12,v15, cylRadius};

#declare edge31 = cylinder{v13,v14, cylRadius};
#declare edge32 = cylinder{v14,v15, cylRadius};


#declare XAxis  = 
object {
     cylinder{<0,0,0> ,<1,0,0> , cylRadius}
     pigment { color Red }
};

#declare YAxis  = 
object {
     cylinder{<0,0,0> ,<0,1,0> , cylRadius}
     pigment { color Green }
};

#declare ZAxis  = 
object {
     cylinder{<0,0,0> ,<0,0,1> , cylRadius}
     pigment { color Blue }
};

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
       object{node9}
       object{node10}
       object{node11}
       object{node12}
       object{node13}
       object{node14}
       object{node15}
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
       object{edge9}
       object{edge10}
       object{edge11}
       object{edge12}

       object{edge13}
       object{edge14}
       object{edge15}
       object{edge16}
       object{edge17}
       object{edge18}
       object{edge19}
       object{edge20}

       object{edge21}
       object{edge22}
       object{edge23}
       object{edge24}
       object{edge25}
       object{edge26}
       object{edge27}
       object{edge28}
       object{edge29}
       object{edge30}
       object{edge31}
       object{edge32}

       texture {T_Silver_1A}
}
  
#declare Scene = 
union {
     object{ XAxis }
     object{ YAxis }
     object{ ZAxis }
     object {HexNodes}
     object {HexEdges}
     rotate<0,360*clock,0>
}
Scene
