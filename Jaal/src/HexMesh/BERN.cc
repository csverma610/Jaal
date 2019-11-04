#include "Mesh.hpp"
#include "tfiblend.hpp"

using namespace Jaal;

void BernHexOps::TestOp314_5()
{
     Mesh *mesh = new Mesh3D;

     JNodeSequence nodes(16);
     for( int i = 0; i < 16; i++) {
          nodes[i] = Vertex::newObject();
          nodes[i]->setID(i);
          mesh->addObject( nodes[i] );
     }
     Point3D p3d;

     p3d[0] = -0.0;
     p3d[1] = -0.0;
     p3d[2] = -0.0;
     nodes[0]->setXYZCoords(p3d);

     // First Plane 
     p3d[0] =  1.0;
     p3d[1] = -0.2;
     p3d[2] =  0.0;
     nodes[1]->setXYZCoords(p3d);
    
     p3d[0] =  1.0;
     p3d[1] =  1.0;
     p3d[2] =  0.0;
     nodes[2]->setXYZCoords(p3d);

     p3d[0] = -0.0;
     p3d[1] =  1.5;
     p3d[2] = -0.0;
     nodes[3]->setXYZCoords(p3d);

     p3d[0] = -1.0;
     p3d[1] =  1.0;
     p3d[2] =  0.0;
     nodes[4]->setXYZCoords(p3d);

     p3d[0] = -1.0;
     p3d[1] = -0.2;
     p3d[2] =  0.0;
     nodes[5]->setXYZCoords(p3d);

     p3d[0] =  0.0;
     p3d[1] = -1.0;
     p3d[2] =  0.0;
     nodes[6]->setXYZCoords(p3d);

     p3d[0] = -0.0;
     p3d[1] = -0.0;
     p3d[2] =  1.0;
     nodes[7]->setXYZCoords(p3d);

     // First Plane 
     p3d[0] =  1.0;
     p3d[1] = -0.2;
     p3d[2] =  1.0;
     nodes[8]->setXYZCoords(p3d);
    
     p3d[0] =  1.0;
     p3d[1] =  1.0;
     p3d[2] =  1.0;
     nodes[9]->setXYZCoords(p3d);

     p3d[0] = -0.0;
     p3d[1] =  1.5;
     p3d[2] =  1.0;
     nodes[10]->setXYZCoords(p3d);

     p3d[0] = -1.0;
     p3d[1] =  1.0;
     p3d[2] =  1.0;
     nodes[11]->setXYZCoords(p3d);

     p3d[0] = -1.0;
     p3d[1] = -0.2;
     p3d[2] =  1.0;
     nodes[12]->setXYZCoords(p3d);

     p3d[0] =  0.0;
     p3d[1] = -1.0;
     p3d[2] =  1.0;
     nodes[13]->setXYZCoords(p3d);

     p3d[0] =  0.0;
     p3d[1] = -0.0;
     p3d[2] =  0.9;
     nodes[14]->setXYZCoords(p3d);

     p3d[0] =  0.0;
     p3d[1] = -0.0;
     p3d[2] =  0.1;
     nodes[15]->setXYZCoords(p3d);


     vector<Vertex*>  hnodes(8);
     hnodes[0] = nodes[0];
     hnodes[1] = nodes[1];
     hnodes[2] = nodes[2];
     hnodes[3] = nodes[3];
     hnodes[4] = nodes[7];
     hnodes[5] = nodes[8];
     hnodes[6] = nodes[9];
     hnodes[7] = nodes[10];
     Cell *hex1 = Hexahedron::newObject();
     hex1->setNodes( hnodes );

     hnodes[0] = nodes[0];
     hnodes[1] = nodes[3];
     hnodes[2] = nodes[4];
     hnodes[3] = nodes[5];
     hnodes[4] = nodes[7];
     hnodes[5] = nodes[10];
     hnodes[6] = nodes[11];
     hnodes[7] = nodes[12];
     Cell *hex2 = Hexahedron::newObject();
     hex2->setNodes( hnodes );

     hnodes[0] = nodes[0];
     hnodes[1] = nodes[5];
     hnodes[2] = nodes[6];
     hnodes[3] = nodes[1];
     hnodes[4] = nodes[7];
     hnodes[5] = nodes[12];
     hnodes[6] = nodes[13];
     hnodes[7] = nodes[8];
     Cell *hex3 = Hexahedron::newObject();
     hex3->setNodes( hnodes );

    mesh->addObject(hex1);
    mesh->addObject(hex2);
    mesh->addObject(hex3);
    mesh->collect_edges();

/*
    mesh->deleteCells();

     hnodes[0] = nodes[9];
     hnodes[1] = nodes[14];
     hnodes[2] = nodes[15];
     hnodes[3] = nodes[2];
     hnodes[4] = nodes[8];
     hnodes[5] = nodes[13];
     hnodes[6] = nodes[6];
     hnodes[7] = nodes[1];
     Cell *hex4 = Hexahedron::newObject();
     hex4->setNodes( hnodes );
     mesh->addObject(hex4);

     hnodes[0] = nodes[11];
     hnodes[1] = nodes[14];
     hnodes[2] = nodes[15];
     hnodes[3] = nodes[4];
     hnodes[4] = nodes[12];
     hnodes[5] = nodes[13];
     hnodes[6] = nodes[6];
     hnodes[7] = nodes[5];
     Cell *hex5 = Hexahedron::newObject();
     hex5->setNodes( hnodes );
     mesh->addObject(hex5);

     hnodes[0] = nodes[11];
     hnodes[1] = nodes[14];
     hnodes[2] = nodes[15];
     hnodes[3] = nodes[4];
     hnodes[4] = nodes[10];
     hnodes[5] = nodes[9];
     hnodes[6] = nodes[2];
     hnodes[7] = nodes[3];
     Cell *hex6 = Hexahedron::newObject();
     hex6->setNodes( hnodes );
     mesh->addObject(hex6);
*/


/*
     Face *face1;
     face1 = Quadrilateral::newObject( nodes[9], nodes[14], nodes[13], nodes[8] );
     mesh->addObject(face1);
     face1 = Quadrilateral::newObject( nodes[11], nodes[14], nodes[13], nodes[12] );
     mesh->addObject(face1);

     face1 = Quadrilateral::newObject( nodes[4], nodes[15], nodes[6], nodes[5] );
     mesh->addObject(face1);
     face1 = Quadrilateral::newObject( nodes[2], nodes[15], nodes[6], nodes[1] );
     mesh->addObject(face1);
*/
     mesh->saveAs( "bern.xml");
}


void BernHexOps::TestOp316_5()
{
     Mesh *mesh = new Mesh3D;

     JNodeSequence nodes(16);
     for( int i = 0; i < 16; i++) {
          nodes[i] = Vertex::newObject();
          nodes[i]->setID(i);
          mesh->addObject( nodes[i] );
     }
     Point3D p3d;

     p3d[0] = -1.0;
     p3d[1] = -1.0;
     p3d[2] =  0.0;
     nodes[0]->setXYZCoords(p3d);

     // First Plane 
     p3d[0] =  1.0;
     p3d[1] = -1.0;
     p3d[2] =  0.0;
     nodes[1]->setXYZCoords(p3d);
    
     p3d[0] =  1.0;
     p3d[1] =  1.0;
     p3d[2] =  0.0;
     nodes[2]->setXYZCoords(p3d);

     p3d[0] = -1.0;
     p3d[1] =  1.0;
     p3d[2] =  0.0;
     nodes[3]->setXYZCoords(p3d);

     p3d[0] = -2.0;
     p3d[1] = -2.0;
     p3d[2] =  1.0;
     nodes[4]->setXYZCoords(p3d);

     p3d[0] =  2.0;
     p3d[1] = -2.0;
     p3d[2] =  1.0;
     nodes[5]->setXYZCoords(p3d);

     p3d[0] =  2.0;
     p3d[1] =  2.0;
     p3d[2] =  1.0;
     nodes[6]->setXYZCoords(p3d);

     p3d[0] = -2.0;
     p3d[1] =  2.0;
     p3d[2] =  1.0;
     nodes[7]->setXYZCoords(p3d);

     p3d[0] = -2.0;
     p3d[1] = -2.0;
     p3d[2] =  2.0;
     nodes[8]->setXYZCoords(p3d);

     p3d[0] =  2.0;
     p3d[1] = -2.0;
     p3d[2] =  2.0;
     nodes[9]->setXYZCoords(p3d);

     p3d[0] =  2.0;
     p3d[1] =  2.0;
     p3d[2] =  2.0;
     nodes[10]->setXYZCoords(p3d);

     p3d[0] = -2.0;
     p3d[1] =  2.0;
     p3d[2] =  2.0;
     nodes[11]->setXYZCoords(p3d);

     p3d[0] = -1.0;
     p3d[1] = -1.0;
     p3d[2] =  3.0;
     nodes[12]->setXYZCoords(p3d);

     // First Plane 
     p3d[0] =  1.0;
     p3d[1] = -1.0;
     p3d[2] =  3.0;
     nodes[13]->setXYZCoords(p3d);
    
     p3d[0] =  1.0;
     p3d[1] =  1.0;
     p3d[2] =  3.0;
     nodes[14]->setXYZCoords(p3d);

     p3d[0] = -1.0;
     p3d[1] =  1.0;
     p3d[2] =  3.0;
     nodes[15]->setXYZCoords(p3d);

     vector<Vertex*>  hnodes(8);
     hnodes[0] = nodes[0];
     hnodes[1] = nodes[1];
     hnodes[2] = nodes[2];
     hnodes[3] = nodes[3];
     hnodes[4] = nodes[4];
     hnodes[5] = nodes[5];
     hnodes[6] = nodes[6];
     hnodes[7] = nodes[7];
     Cell *hex1 = Hexahedron::newObject();
     hex1->setNodes( hnodes );

     hnodes[0] = nodes[4];
     hnodes[1] = nodes[5];
     hnodes[2] = nodes[6];
     hnodes[3] = nodes[7];
     hnodes[4] = nodes[8];
     hnodes[5] = nodes[9];
     hnodes[6] = nodes[10];
     hnodes[7] = nodes[11];
     Cell *hex2 = Hexahedron::newObject();
     hex2->setNodes( hnodes );

     hnodes[0] = nodes[8];
     hnodes[1] = nodes[9];
     hnodes[2] = nodes[10];
     hnodes[3] = nodes[11];
     hnodes[4] = nodes[12];
     hnodes[5] = nodes[13];
     hnodes[6] = nodes[14];
     hnodes[7] = nodes[15];
     Cell *hex3 = Hexahedron::newObject();
     hex3->setNodes( hnodes );

     mesh->addObject(hex1);
     mesh->addObject(hex2);
     mesh->addObject(hex3);
/*
     face1 = Quadrilateral::newObject( nodes[0], nodes[1], nodes[13], nodes[12] );
     mesh->addObject(face1);
     face1 = Quadrilateral::newObject( nodes[3], nodes[2], nodes[14], nodes[15] );
     mesh->addObject(face1);
     face1 = Quadrilateral::newObject( nodes[8], nodes[9], nodes[13], nodes[12] );
     mesh->addObject(face1);
     face1 = Quadrilateral::newObject( nodes[11], nodes[10], nodes[14], nodes[15] );
     mesh->addObject(face1);
*/

     mesh->collect_edges();

     mesh->saveAs( "bern.xml");
}


void BernHexOps::TestOp416_4()
{
     Mesh *mesh = new Mesh3D;

     JNodeSequence nodes(16);
     for( int i = 0; i < 16; i++) {
          nodes[i] = Vertex::newObject();
          nodes[i]->setID(i);
          mesh->addObject( nodes[i] );
     }
     Point3D p3d;

     p3d[0] = -1.0;
     p3d[1] = -0.0;
     p3d[2] = -0.1;
     nodes[0]->setXYZCoords(p3d);

     // First Plane 
     p3d[0] =  1.0;
     p3d[1] =  0.0;
     p3d[2] = -0.1;
     nodes[1]->setXYZCoords(p3d);
    
     p3d[0] =  1.5;
     p3d[1] =  1.0;
     p3d[2] =  0.0;
     nodes[2]->setXYZCoords(p3d);

     p3d[0] = -1.5;
     p3d[1] =  1.0;
     p3d[2] = 0.0;
     nodes[3]->setXYZCoords(p3d);

     p3d[0] = -1.5;
     p3d[1] = -1.0;
     p3d[2] =  0.0;
     nodes[4]->setXYZCoords(p3d);

     p3d[0] =  1.5;
     p3d[1] = -1.0;
     p3d[2] = -0.0;
     nodes[5]->setXYZCoords(p3d);

     p3d[0] =  2.0;
     p3d[1] =  0.0;
     p3d[2] = -0.0;
     nodes[6]->setXYZCoords(p3d);

     p3d[0] = -2.0;
     p3d[1] =  0.0;
     p3d[2] = -0.0;
     nodes[7]->setXYZCoords(p3d);

     // First Plane 
     p3d[0] = -1.0;
     p3d[1] =  0.0;
     p3d[2] =  1.1;
     nodes[8]->setXYZCoords(p3d);

     p3d[0] =  1.0;
     p3d[1] = -0.0;
     p3d[2] =  1.1;
     nodes[9]->setXYZCoords(p3d);
    
     p3d[0] =  1.5;
     p3d[1] =  1.0;
     p3d[2] =  1.0;
     nodes[10]->setXYZCoords(p3d);

     p3d[0] = -1.5;
     p3d[1] =  1.0;
     p3d[2] =  1.0;
     nodes[11]->setXYZCoords(p3d);

     p3d[0] = -1.5;
     p3d[1] = -1.0;
     p3d[2] =  1.0;
     nodes[12]->setXYZCoords(p3d);

     p3d[0] =  1.5;
     p3d[1] = -1.0;
     p3d[2] =  1.0;
     nodes[13]->setXYZCoords(p3d);

     p3d[0] =  2.0;
     p3d[1] =  0.0;
     p3d[2] =  1.0;
     nodes[14]->setXYZCoords(p3d);

     p3d[0] = -2.0;
     p3d[1] =  0.0;
     p3d[2] =  1.0;
     nodes[15]->setXYZCoords(p3d);


     vector<Vertex*>  hnodes(8);
     hnodes[0] = nodes[0];
     hnodes[1] = nodes[1];
     hnodes[2] = nodes[2];
     hnodes[3] = nodes[3];
     hnodes[4] = nodes[8];
     hnodes[5] = nodes[9];
     hnodes[6] = nodes[10];
     hnodes[7] = nodes[11];
     Cell *hex1 = Hexahedron::newObject();
     hex1->setNodes( hnodes );

     hnodes[0] = nodes[1];
     hnodes[1] = nodes[0];
     hnodes[2] = nodes[4];
     hnodes[3] = nodes[5];
     hnodes[4] = nodes[9];
     hnodes[5] = nodes[8];
     hnodes[6] = nodes[12];
     hnodes[7] = nodes[13];
     Cell *hex2 = Hexahedron::newObject();
     hex2->setNodes( hnodes );

     hnodes[0] = nodes[1];
     hnodes[1] = nodes[5];
     hnodes[2] = nodes[6];
     hnodes[3] = nodes[2];
     hnodes[4] = nodes[9];
     hnodes[5] = nodes[13];
     hnodes[6] = nodes[14];
     hnodes[7] = nodes[10];
     Cell *hex3 = Hexahedron::newObject();
     hex3->setNodes( hnodes );

     hnodes[0] = nodes[3];
     hnodes[1] = nodes[7];
     hnodes[2] = nodes[4];
     hnodes[3] = nodes[0];
     hnodes[4] = nodes[11];
     hnodes[5] = nodes[15];
     hnodes[6] = nodes[12];
     hnodes[7] = nodes[8];
     Cell *hex4 = Hexahedron::newObject();
     hex4->setNodes( hnodes );

     mesh->addObject(hex1);
     mesh->addObject(hex2);
     mesh->addObject(hex3);
     mesh->addObject(hex4);
/*
     face1 = Quadrilateral::newObject( nodes[0], nodes[1], nodes[13], nodes[12] );
     mesh->addObject(face1);
     face1 = Quadrilateral::newObject( nodes[3], nodes[2], nodes[14], nodes[15] );
     mesh->addObject(face1);
     face1 = Quadrilateral::newObject( nodes[8], nodes[9], nodes[13], nodes[12] );
     mesh->addObject(face1);
     face1 = Quadrilateral::newObject( nodes[11], nodes[10], nodes[14], nodes[15] );
     mesh->addObject(face1);
*/

     mesh->collect_edges();

     mesh->saveAs( "bern.xml");
}

void BernHexOps::TestOp2_6()
{
     Mesh *mesh = new Mesh3D;

     JNodeSequence nodes(16);
     for( int i = 0; i < 16; i++) {
          nodes[i] = Vertex::newObject();
          nodes[i]->setID(i);
          mesh->addObject( nodes[i] );
     }
     Point3D p3d;

     // First Plane 
     p3d[0] = -1.0;
     p3d[1] = -1.0;
     p3d[2] = -2.0;
     nodes[0]->setXYZCoords(p3d);
    
     p3d[0] =  1.0;
     p3d[1] = -1.0;
     p3d[2] = -2.0;
     nodes[1]->setXYZCoords(p3d);

     p3d[0] =  1.0;
     p3d[1] =  1.0;
     p3d[2] = -2.0;
     nodes[2]->setXYZCoords(p3d);

     p3d[0] = -1.0;
     p3d[1] =  1.0;
     p3d[2] = -2.0;
     nodes[3]->setXYZCoords(p3d);

     // Second Plane 
     p3d[0] = -2.0;
     p3d[1] = -2.0;
     p3d[2] = -0.0;
     nodes[4]->setXYZCoords(p3d);
    
     p3d[0] =  2.0;
     p3d[1] = -2.0;
     p3d[2] = -0.0;
     nodes[5]->setXYZCoords(p3d);

     p3d[0] =  2.0;
     p3d[1] =  2.0;
     p3d[2] = -0.0;
     nodes[6]->setXYZCoords(p3d);

     p3d[0] = -2.0;
     p3d[1] =  2.0;
     p3d[2] =  0.0;
     nodes[7]->setXYZCoords(p3d);

     // Third Plane
     p3d[0] = -1.0;
     p3d[1] = -1.0;
     p3d[2] =  2.0;
     nodes[8]->setXYZCoords(p3d);
    
     p3d[0] =  1.0;
     p3d[1] = -1.0;
     p3d[2] =  2.0;
     nodes[9]->setXYZCoords(p3d);

     p3d[0] =  1.0;
     p3d[1] =  1.0;
     p3d[2] =  2.0;
     nodes[10]->setXYZCoords(p3d);

     p3d[0] = -1.0;
     p3d[1] =  1.0;
     p3d[2] =  2.0;
     nodes[11]->setXYZCoords(p3d);

     // MidPlane 
     p3d[0] = -1.0;
     p3d[1] = -1.0;
     p3d[2] = -0.0;
     nodes[12]->setXYZCoords(p3d);
    
     p3d[0] =  1.0;
     p3d[1] = -1.0;
     p3d[2] = -0.0;
     nodes[13]->setXYZCoords(p3d);

     p3d[0] =  1.0;
     p3d[1] =  1.0;
     p3d[2] = -0.0;
     nodes[14]->setXYZCoords(p3d);

     p3d[0] = -1.0;
     p3d[1] =  1.0;
     p3d[2] =  0.0;
     nodes[15]->setXYZCoords(p3d);

     Face *face1 = Quadrilateral::newObject( nodes[12], nodes[13], nodes[14], nodes[15] );
     mesh->addObject( face1 );

     vector<Vertex*>  hnodes(8);
     hnodes[0] = nodes[0];
     hnodes[1] = nodes[1];
     hnodes[2] = nodes[2];
     hnodes[3] = nodes[3];
     hnodes[4] = nodes[4];
     hnodes[5] = nodes[5];
     hnodes[6] = nodes[6];
     hnodes[7] = nodes[7];
     Cell *hex1 = Hexahedron::newObject();
     hex1->setNodes( hnodes );
     hnodes[0] = nodes[4];
     hnodes[1] = nodes[5];
     hnodes[2] = nodes[6];
     hnodes[3] = nodes[7];
     hnodes[4] = nodes[8];
     hnodes[5] = nodes[9];
     hnodes[6] = nodes[10];
     hnodes[7] = nodes[11];
     Cell *hex2 = Hexahedron::newObject();
     hex2->setNodes( hnodes );

     mesh->addObject(hex1);
     mesh->addObject(hex2);
     mesh->saveAs( "bern1.xml");

     face1 = Quadrilateral::newObject( nodes[0], nodes[1], nodes[13], nodes[12] );
     mesh->addObject(face1);

     face1 = Quadrilateral::newObject( nodes[3], nodes[2], nodes[14], nodes[15] );
     mesh->addObject(face1);

     face1 = Quadrilateral::newObject( nodes[8], nodes[9], nodes[13], nodes[12] );
     mesh->addObject(face1);

     face1 = Quadrilateral::newObject( nodes[11], nodes[10], nodes[14], nodes[15] );
     mesh->addObject(face1);
 
    mesh->deleteCells();
    Hexahedron *h1 = Hexahedron::newObject();
    h1->setNodes( nodes[0], nodes[1], nodes[2], nodes[3], nodes[12], nodes[13], nodes[14], nodes[16]);
                

}

///////////////////////////////////////////////////////////////////////////////////

void BernHexOps::Op1_7( const Hexahedron *hex, JNodeSequence &newnodes, JCellSequence &newcells)
{
     newcells.resize(7);
     double xc[8], yc[8], zc[8];

     for( int i = 0; i < 8; i++) {
         const Point3D &p3d = hex->getNodeAt(i)->getXYZCoords();
         xc[i] = p3d[0];
         yc[i] = p3d[1];
         zc[i] = p3d[2];
     }

     newnodes.resize(8);

     Vertex *vtx;
     Point3D p3d;

     vtx = Vertex::newObject();
     p3d[0] = TFI::trilinear_interpolation( -0.5, -0.5, -0.5, xc);
     p3d[1] = TFI::trilinear_interpolation( -0.5, -0.5, -0.5, yc);
     p3d[2] = TFI::trilinear_interpolation( -0.5, -0.5, -0.5, zc);
     vtx->setXYZCoords(p3d);
     newnodes[0]  = vtx;
// CSV  

     vtx = Vertex::newObject();
     p3d[0] = TFI::trilinear_interpolation( 0.5, -0.5, -0.5, xc);
     p3d[1] = TFI::trilinear_interpolation( 0.5, -0.5, -0.5, yc);
     p3d[2] = TFI::trilinear_interpolation( 0.5, -0.5, -0.5, zc);
     vtx->setXYZCoords(p3d);
     newnodes[1]  = vtx;

     vtx = Vertex::newObject();
     p3d[0] = TFI::trilinear_interpolation( 0.5, 0.5, -0.5, xc);
     p3d[1] = TFI::trilinear_interpolation( 0.5, 0.5, -0.5, yc);
     p3d[2] = TFI::trilinear_interpolation( 0.5, 0.5, -0.5, zc);
     vtx->setXYZCoords(p3d);
     newnodes[2]  = vtx;

     vtx = Vertex::newObject();
     p3d[0] = TFI::trilinear_interpolation( -0.5, 0.5, -0.5, xc);
     p3d[1] = TFI::trilinear_interpolation( -0.5, 0.5, -0.5, yc);
     p3d[2] = TFI::trilinear_interpolation( -0.5, 0.5, -0.5, zc);
     vtx->setXYZCoords(p3d);
     newnodes[3]  = vtx;

     vtx = Vertex::newObject();
     p3d[0] = TFI::trilinear_interpolation( -0.5, -0.5, 0.5, xc);
     p3d[1] = TFI::trilinear_interpolation( -0.5, -0.5, 0.5, yc);
     p3d[2] = TFI::trilinear_interpolation( -0.5, -0.5, 0.5, zc);
     vtx->setXYZCoords(p3d);
     newnodes[4]  = vtx;

     vtx = Vertex::newObject();
     p3d[0] = TFI::trilinear_interpolation( 0.5, -0.5, 0.5, xc);
     p3d[1] = TFI::trilinear_interpolation( 0.5, -0.5, 0.5, yc);
     p3d[2] = TFI::trilinear_interpolation( 0.5, -0.5, 0.5, zc);
     vtx->setXYZCoords(p3d);
     newnodes[5]  = vtx;

     vtx = Vertex::newObject();
     p3d[0] = TFI::trilinear_interpolation( 0.5, 0.5, 0.5, xc);
     p3d[1] = TFI::trilinear_interpolation( 0.5, 0.5, 0.5, yc);
     p3d[2] = TFI::trilinear_interpolation( 0.5, 0.5, 0.5, zc);
     vtx->setXYZCoords(p3d);
     newnodes[6]  = vtx;

     vtx = Vertex::newObject();
     p3d[0] = TFI::trilinear_interpolation( -0.5, 0.5, 0.5, xc);
     p3d[1] = TFI::trilinear_interpolation( -0.5, 0.5, 0.5, yc);
     p3d[2] = TFI::trilinear_interpolation( -0.5, 0.5, 0.5, zc);
     vtx->setXYZCoords(p3d);
     newnodes[7]  = vtx;

     Hexahedron *hex0 = Hexahedron::newObject();
     hex0->setNodes(newnodes);

     JNodeSequence hexnodes(8);

     Hexahedron *hex1 = Hexahedron::newObject();
     hexnodes[0] =  hex->getNodeAt(0);
     hexnodes[1] =  hex->getNodeAt(4);
     hexnodes[2] =  hex->getNodeAt(7);
     hexnodes[3] =  hex->getNodeAt(3);
     hexnodes[4] =  newnodes[0];
     hexnodes[5] =  newnodes[4];
     hexnodes[6] =  newnodes[7];
     hexnodes[7] =  newnodes[3];
     hex1->setNodes( hexnodes );

     Hexahedron *hex2 = Hexahedron::newObject();
     hexnodes[0] =  hex->getNodeAt(1);
     hexnodes[1] =  hex->getNodeAt(5);
     hexnodes[2] =  hex->getNodeAt(6);
     hexnodes[3] =  hex->getNodeAt(2);
     hexnodes[4] =  newnodes[1];
     hexnodes[5] =  newnodes[5];
     hexnodes[6] =  newnodes[6];
     hexnodes[7] =  newnodes[2];
     hex2->setNodes( hexnodes );

     Hexahedron *hex3 = Hexahedron::newObject();
     hexnodes[0] =  hex->getNodeAt(0);
     hexnodes[1] =  hex->getNodeAt(1);
     hexnodes[2] =  hex->getNodeAt(5);
     hexnodes[3] =  hex->getNodeAt(4);
     hexnodes[4] =  newnodes[0];
     hexnodes[5] =  newnodes[1];
     hexnodes[6] =  newnodes[5];
     hexnodes[7] =  newnodes[4];
     hex3->setNodes( hexnodes );


     Hexahedron *hex4 = Hexahedron::newObject();
     hexnodes[0] =  hex->getNodeAt(3);
     hexnodes[1] =  hex->getNodeAt(2);
     hexnodes[2] =  hex->getNodeAt(6);
     hexnodes[3] =  hex->getNodeAt(7);
     hexnodes[4] =  newnodes[3];
     hexnodes[5] =  newnodes[2];
     hexnodes[6] =  newnodes[6];
     hexnodes[7] =  newnodes[7];
     hex4->setNodes( hexnodes );

     Hexahedron *hex5 = Hexahedron::newObject();
     hexnodes[0] =  hex->getNodeAt(0);
     hexnodes[1] =  hex->getNodeAt(1);
     hexnodes[2] =  hex->getNodeAt(2);
     hexnodes[3] =  hex->getNodeAt(3);
     hexnodes[4] =  newnodes[0];
     hexnodes[5] =  newnodes[1];
     hexnodes[6] =  newnodes[2];
     hexnodes[7] =  newnodes[3];
     hex5->setNodes( hexnodes );

     Hexahedron *hex6 = Hexahedron::newObject();
     hexnodes[0] =  hex->getNodeAt(4);
     hexnodes[1] =  hex->getNodeAt(5);
     hexnodes[2] =  hex->getNodeAt(6);
     hexnodes[3] =  hex->getNodeAt(7);
     hexnodes[4] =  newnodes[4];
     hexnodes[5] =  newnodes[5];
     hexnodes[6] =  newnodes[6];
     hexnodes[7] =  newnodes[7];
     hex6->setNodes( hexnodes );

     newcells.resize(7);
     newcells[0] = hex0;
     newcells[1] = hex1;
     newcells[2] = hex2;
     newcells[3] = hex3;
     newcells[4] = hex4;
     newcells[5] = hex5;
     newcells[6] = hex6;
}
////////////////////////////////////////////////////////////////////////////////

void BernHexOps::Op2_6( const Hexahedron *hex1, const Hexahedron *hex2, JNodeSequence &newnodes, JCellSequence &newcells)
{
     newcells.resize(6);
     JNodeSequence nodes(16), hexnodes(8);

     newnodes.resize(4);
     newnodes[0] = nodes[12];
     newnodes[1] = nodes[13];
     newnodes[2] = nodes[14];
     newnodes[3] = nodes[15];

      // First Cell...
     hexnodes[0] = nodes[0];
     hexnodes[1] = nodes[1];
     hexnodes[2] = nodes[2];
     hexnodes[3] = nodes[3];
     hexnodes[4] = nodes[12];
     hexnodes[5] = nodes[13];
     hexnodes[6] = nodes[14];
     hexnodes[7] = nodes[15];
     newcells[0] = Hexahedron::newObject();
     newcells[0]->setNodes( hexnodes );

     hexnodes[0] = nodes[12];
     hexnodes[1] = nodes[13];
     hexnodes[2] = nodes[14];
     hexnodes[3] = nodes[15];
     hexnodes[4] = nodes[8];
     hexnodes[5] = nodes[9];
     hexnodes[6] = nodes[10];
     hexnodes[7] = nodes[11];
     newcells[1] = Hexahedron::newObject();
     newcells[1]->setNodes( hexnodes );

     hexnodes[0] = nodes[1];
     hexnodes[1] = nodes[13];
     hexnodes[2] = nodes[9];
     hexnodes[3] = nodes[5];
     hexnodes[4] = nodes[2];
     hexnodes[5] = nodes[14];
     hexnodes[6] = nodes[10];
     hexnodes[7] = nodes[6];
     newcells[2] = Hexahedron::newObject();
     newcells[2]->setNodes( hexnodes );

     hexnodes[0] = nodes[3];
     hexnodes[1] = nodes[15];
     hexnodes[2] = nodes[11];
     hexnodes[3] = nodes[7];
     hexnodes[4] = nodes[0];
     hexnodes[5] = nodes[12];
     hexnodes[6] = nodes[8];
     hexnodes[7] = nodes[4];
     newcells[3] = Hexahedron::newObject();
     newcells[3]->setNodes( hexnodes );

     hexnodes[0] = nodes[0];
     hexnodes[1] = nodes[12];
     hexnodes[2] = nodes[8];
     hexnodes[3] = nodes[4];
     hexnodes[4] = nodes[1];
     hexnodes[5] = nodes[13];
     hexnodes[6] = nodes[9];
     hexnodes[7] = nodes[5];
     newcells[4] = Hexahedron::newObject();
     newcells[4]->setNodes( hexnodes );
     
     hexnodes[0] = nodes[3];
     hexnodes[1] = nodes[15];
     hexnodes[2] = nodes[11];
     hexnodes[3] = nodes[7];
     hexnodes[4] = nodes[2];
     hexnodes[5] = nodes[14];
     hexnodes[6] = nodes[10];
     hexnodes[7] = nodes[6];
     newcells[5] = Hexahedron::newObject();
     newcells[5]->setNodes( hexnodes );
}

////////////////////////////////////////////////////////////////////////////////
void BernHexOps::Op314_5( const Hexahedron *hex1, const Hexahedron *hex2, const Hexahedron *hex3, JCellSequence &newCells)
{

    newCells.resize(5);

     /*
     hexnodes[0] = nodes[9];
     hexnodes[1] = nodes[14];
     hexnodes[2] = nodes[15];
     hexnodes[3] = nodes[2];
     hexnodes[4] = nodes[8];
     hexnodes[5] = nodes[13];
     hexnodes[6] = nodes[6];
     hexnodes[7] = nodes[1];
     newcells[0] = Hexahedron::newObject();
     newcells[0]->setNodes( hexnodes );

     hexnodes[0] = nodes[11];
     hexnodes[1] = nodes[14];
     hexnodes[2] = nodes[15];
     hexnodes[3] = nodes[4];
     hexnodes[4] = nodes[12];
     hexnodes[5] = nodes[13];
     hexnodes[6] = nodes[6];
     hexnodes[7] = nodes[5];
     newcells[1] = Hexahedron::newObject();
     newcells[1]->setNodes( hexnodes );

     hexnodes[0] = nodes[11];
     hexnodes[1] = nodes[14];
     hexnodes[2] = nodes[15];
     hexnodes[3] = nodes[4];
     hexnodes[4] = nodes[10];
     hexnodes[5] = nodes[9];
     hexnodes[6] = nodes[2];
     hexnodes[7] = nodes[3];
     newcells[3] = Hexahedron::newObject();
     newcells[3]->setNodes( hexnodes );

     hexnodes[0] = nodes[10];
     hexnodes[1] = nodes[9];
     hexnodes[2] = nodes[8];
     hexnodes[3] = nodes[7];
     hexnodes[4] = nodes[11];
     hexnodes[5] = nodes[14];
     hexnodes[6] = nodes[13];
     hexnodes[7] = nodes[12];
     newcells[4] = Hexahedron::newObject();
     newcells[4]->setNodes( hexnodes );

     hexnodes[0] = nodes[3];
     hexnodes[1] = nodes[2];
     hexnodes[2] = nodes[1];
     hexnodes[3] = nodes[0];
     hexnodes[4] = nodes[4];
     hexnodes[5] = nodes[15];
     hexnodes[6] = nodes[6];
     hexnodes[7] = nodes[5];
     newcells[5] = Hexahedron::newObject();
     newcells[5]->setNodes( hexnodes );
    */
}
///////////////////////////////////////////////////////////////////////////////

void BernHexOps:: Op316_5( const Hexahedron *hex1, const Hexahedron *hex2, const Hexahedron *hex3, JCellSequence &newCells)
{
/*
    Face tCommFace, bCommFace;  // Top and Bottom Common faces.
    Face tFace, bFace;          // Top and Bottom faces.

    int stat = Hexahedron::shared_entitity( tHex, mHex, tCommFace);
    if( stat == 0) return;

    int stat = Hexahedron::shared_entitity( bHex, mHex, bCommFace);
    if( stat == 0) return;

    tHex->getOppositeFace( tCommFace, tFace );   // Face opposite to the Top  common face ...
    bHex->getOppositeFace( bCommFace, bFace );   // Face opposite to the Bottom common face ...
*/

/*
   newCells.resize(5);
   JNodeSequence hexnodes(8);

   hexnodes[0] = 0;
   hexnodes[1] = 1;
   hexnodes[2] = 2;
   hexnodes[3] = 3;
   hexnodes[4] = 12;
   hexnodes[5] = 13;
   hexnodes[6] = 14;
   hexnodes[7] = 15;
   newCells[0] = Hexahedron::newObject(hexnodes);

   hexnodes[0] = 1;
   hexnodes[1] = 5;
   hexnodes[2] = 6;
   hexnodes[3] = 2;
   hexnodes[4] = 13;
   hexnodes[5] =  9;
   hexnodes[6] = 10;
   hexnodes[7] = 14;
   newCells[1] = Hexahedron::newObject(hexnodes);
   
   hexnodes[0] = 0;
   hexnodes[1] = 3;
   hexnodes[2] = 7;
   hexnodes[3] = 4;
   hexnodes[4] = 12;
   hexnodes[5] = 15;
   hexnodes[6] = 11;
   hexnodes[7] = 8;
   newCells[2] = Hexahedron::newObject(hexnodes);

   hexnodes[0] = 3;
   hexnodes[1] = 2;
   hexnodes[2] = 6;
   hexnodes[3] = 7;
   hexnodes[4] = 15;
   hexnodes[5] = 14;
   hexnodes[6] = 10;
   hexnodes[7] = 11;
   newCells[3] = Hexahedron::newObject(hexnodes);

   hexnodes[0] = 0;
   hexnodes[1] = 4;
   hexnodes[2] = 5;
   hexnodes[3] = 1;
   hexnodes[4] = 12;
   hexnodes[5] = 8;
   hexnodes[6] = 9;
   hexnodes[7] = 13;
   newCells[4] = Hexahedron::newObject(hexnodes);
*/
   
}

///////////////////////////////////////////////////////////////////////////////////
void BernHexOps:: Op416_4( const Hexahedron *hex1, const Hexahedron *hex2, 
                        const Hexahedron *hex3, const Hexahedron *hex4,
                        JCellSequence &newcells)
{
    JNodeSequence nodes(16), hexnodes(8);

    // Let the common face is (0,1,9,8);

    newcells.resize(4);

     hexnodes[0] = nodes[2];
     hexnodes[1] = nodes[10];
     hexnodes[2] = nodes[14];
     hexnodes[3] = nodes[6];
     hexnodes[4] = nodes[3];
     hexnodes[5] = nodes[11];
     hexnodes[6] = nodes[15];
     hexnodes[7] = nodes[7];
     newcells[0] = Hexahedron::newObject();
     newcells[0]->setNodes( hexnodes );

     hexnodes[0] = nodes[6];
     hexnodes[1] = nodes[14];
     hexnodes[2] = nodes[13];
     hexnodes[3] = nodes[5];
     hexnodes[4] = nodes[7];
     hexnodes[5] = nodes[15];
     hexnodes[6] = nodes[12];
     hexnodes[7] = nodes[4];
     newcells[1] = Hexahedron::newObject();
     newcells[1]->setNodes( hexnodes );

     hexnodes[0] = nodes[9];
     hexnodes[1] = nodes[13];
     hexnodes[2] = nodes[14];
     hexnodes[3] = nodes[10];
     hexnodes[4] = nodes[8];
     hexnodes[5] = nodes[12];
     hexnodes[6] = nodes[15];
     hexnodes[7] = nodes[11];
     newcells[2] = Hexahedron::newObject();
     newcells[2]->setNodes( hexnodes );

     hexnodes[0] = nodes[1];
     hexnodes[1] = nodes[2];
     hexnodes[2] = nodes[6];
     hexnodes[3] = nodes[5];
     hexnodes[4] = nodes[0];
     hexnodes[5] = nodes[3];
     hexnodes[6] = nodes[7];
     hexnodes[7] = nodes[4];
     newcells[3] = Hexahedron::newObject();
     newcells[3]->setNodes( hexnodes );
}

