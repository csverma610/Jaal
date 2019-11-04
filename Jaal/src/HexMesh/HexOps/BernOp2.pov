#include "colors.inc"
#include "woods.inc"
#include "stones.inc"
#include "metals.inc"
#include "golds.inc"
#include "glass.inc"

camera {
  sky <0,0,1>          
  location  <0,  2, -3>  
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

#declare L  = 1.0;

#declare v0  = < -L,0,0 >;
#declare v1  = <  0, -0.5*L, -0.5*L >;
#declare v2  = <  0,  1.5*L, -0.5*L >;
#declare v3  = < -L,L,0 >;

#declare v4  = < -L,0,L >;
#declare v5  = <  0,-0.5*L, 1.5*L >;
#declare v6  = <  0, 1.5*L, 1.5*L >;
#declare v7  = < -L,L,L >;

#declare v8  = <  L, 0, 0>;
#declare v9  = <  L, L, 0>;
#declare v10 = <  L, 0, L>;
#declare v11 = <  L, L, L>;

#declare v12  = <  0,0,0 >;
#declare v13  = <  0,L,0 >;
#declare v14  = <  0,0,L >;
#declare v15  = <  0,L,L >;

//Create a "floor"
plane {
  <0,1,0>,-1.1
  texture { T_Stone5 }
}

#declare cylRadius = 0.02;
#declare sphRadius = 0.05;

#declare node0 = sphere{v0, sphRadius};
#declare node1 = sphere{v1, sphRadius};
#declare node2 = sphere{v2, sphRadius};
#declare node3 = sphere{v3, sphRadius};
#declare node4 = sphere{v4, sphRadius};
#declare node5 = sphere{v5, sphRadius};
#declare node6 = sphere{v6, sphRadius};
#declare node7 = sphere{v7, sphRadius};

#declare node8  = sphere{v8,  sphRadius};
#declare node9  = sphere{v9,  sphRadius};
#declare node10 = sphere{v10, sphRadius};
#declare node11 = sphere{v11, sphRadius};

#declare node12 = sphere{v12, sphRadius};
#declare node13 = sphere{v13, sphRadius};
#declare node14 = sphere{v14, sphRadius};
#declare node15 = sphere{v15, sphRadius};

#declare edge_v0v1  = cylinder{v0,v1, cylRadius};
#declare edge_v0v3  = cylinder{v0,v3, cylRadius};
#declare edge_v0v4  = cylinder{v0,v4, cylRadius};

#declare edge_v1v2  = cylinder{v1,v2, cylRadius};
#declare edge_v1v5  = cylinder{v1,v5, cylRadius};
#declare edge_v2v3  = cylinder{v2,v3, cylRadius};

#declare edge_v2v6  = cylinder{v2,v6, cylRadius};
#declare edge_v3v7  = cylinder{v3,v7, cylRadius};
#declare edge_v4v5  = cylinder{v4,v5, cylRadius};
#declare edge_v4v7  = cylinder{v4,v7, cylRadius};
#declare edge_v5v6  = cylinder{v5,v6, cylRadius};
#declare edge_v6v7  = cylinder{v6,v7, cylRadius};

#declare edge_v1v8  = cylinder{v1,v8, cylRadius};
#declare edge_v2v9  = cylinder{v2,v9, cylRadius};
#declare edge_v5v10 = cylinder{v5,v10, cylRadius};
#declare edge_v6v11 = cylinder{v6,v11, cylRadius};

#declare edge_v8v9   = cylinder{v8,v9, cylRadius};
#declare edge_v8v10  = cylinder{v8,v10, cylRadius};
#declare edge_v9v11  = cylinder{v9,v11, cylRadius};
#declare edge_v10v11 = cylinder{v10,v11, cylRadius};

#declare edge_v0v12 = cylinder{v0,v12, cylRadius};
#declare edge_v3v13 = cylinder{v3,v13, cylRadius};
#declare edge_v4v14 = cylinder{v4,v14, cylRadius};
#declare edge_v7v15 = cylinder{v7,v15, cylRadius};

#declare edge_v8v12  = cylinder{v8,v12, cylRadius};
#declare edge_v9v13  = cylinder{v9,v13, cylRadius};
#declare edge_v10v14 = cylinder{v10,v14, cylRadius};
#declare edge_v11v15 = cylinder{v11,v15, cylRadius};

#declare edge_v12v14  = cylinder{v12,v14, cylRadius};
#declare edge_v12v13  = cylinder{v12,v13, cylRadius};
#declare edge_v14v15  = cylinder{v14,v15, cylRadius};
#declare edge_v13v15  = cylinder{v13,v15, cylRadius};

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

#declare Hex1Nodes = 
union { 
       object{node0}
       object{node1}
       object{node2}
       object{node3}
       object{node4}
       object{node5}
       object{node6}
       object{node7}
       texture {T_Gold_1A}
}

#declare Hex1Edges = 
union { 
       object{edge_v0v1}
       object{edge_v1v2}
       object{edge_v2v3}
       object{edge_v0v3}

       object{edge_v4v5}
       object{edge_v5v6}
       object{edge_v6v7}
       object{edge_v4v7}

       object{edge_v0v4}
       object{edge_v3v7}
       object{edge_v2v6}
       object{edge_v1v5}

       texture {T_Silver_1A}
}

#declare Hex2Nodes = 
union { 
       object{node8}
       object{node9}
       object{node10}
       object{node11}
       texture {T_Gold_1A}
}

#declare Hex2Edges = 
union { 
       object{edge_v1v8}
       object{edge_v8v9}
       object{edge_v2v9}

       object{edge_v5v10}
       object{edge_v10v11}
       object{edge_v6v11}

       object{edge_v8v10}
       object{edge_v9v11}

       texture {T_Silver_1A}
}

#declare BridgeNodes = 
union { 
       object{node12}
       object{node13}
       object{node14}
       object{node15}
       texture {T_Gold_1A}
}

#declare BridgeEdges = 
union { 
       object{edge_v0v12}
       object{edge_v3v13}
       object{edge_v4v14}
       object{edge_v7v15}

       object{edge_v8v12}
       object{edge_v9v13}
       object{edge_v10v14}
       object{edge_v11v15}

       object{edge_v12v14}
       object{edge_v12v13}
       object{edge_v14v15}
       object{edge_v13v15}

       texture {T_Silver_1A}
}

#declare Face_v0_v4_v7_v3 =  
union {
      object{ polygon{4, v0, v4, v7, v3} }
      pigment{ color Red }
      finish { phong 1.0 }
}

#declare Face_v3_v7_v15_v13 =  
union {
      object{ polygon{4, v3, v7, v15, v13} }
      pigment{ color Red }
      finish { phong 1.0 }
}

#declare Face_v12_v13_v15_v14 =  
union {
      object{ polygon{4, v12, v13, v15, v14} }
      pigment{ color Red }
      finish { phong 1.0 }
}

#declare Face4 =  
union {
      object{ polygon{4, v0, v12, v14, v4} }
      pigment{ color Red }
      finish { phong 1.0 }
}

#declare Face5 =  
union {
      object{ polygon{4, v4, v14, v15, v7} }
      pigment{ color Red }
      finish { phong 1.0 }
}

#declare Face6 =  
union {
      object{ polygon{4, v0, v12, v13, v3} }
      pigment{ color Red }
      finish { phong 1.0 }
}

#declare Hex3Faces = 
union {
     object { Face_v0_v4_v7_v3 }
     object { Face_v12_v13_v15_v14 }
     object { Face_v3_v7_v15_v13 }
}

#declare Scene = 
union {
     object {Hex1Nodes}
     object {Hex2Nodes}
     object {Hex1Edges}
     object {Hex2Edges}
     object {BridgeNodes}
     object {BridgeEdges}
     object {Hex3Faces}
     rotate<0,60,0>
}

Scene
