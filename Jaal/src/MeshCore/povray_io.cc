#include "MeshExporter.hpp"
#include "MeshGeometry.hpp"
#include "MeshTopology.hpp"

void JMeshPOVExporter :: writeheader( const JMeshPtr &mesh, ofstream &ofile )
{
     ofile << "#include \"colors.inc\""   << endl;
     ofile << "#include \"textures.inc\"" << endl;
     ofile << "#include \"metals.inc\""   << endl;
     ofile << "#include \"woods.inc\""    << endl;
     ofile << "#include \"skies.inc\""    << endl;
     ofile << "#include \"stones.inc\""   << endl << endl;

     ofile << "global_settings" << endl;
     ofile << "{" << endl;
     ofile << "   ambient_light color rgb <1.0, 1.0, 1.0>" << endl;
     ofile << "   assumed_gamma 2" << endl;
     ofile << "#declare RAD= 0;" << endl;

     ofile << "#if(RAD)" << endl;
     ofile << "   radiosity" << endl;
     ofile << "   {" << endl;
     ofile << "    brightness 2.0" << endl;
     ofile << "    count 10" << endl;
     ofile << "    error_bound 0.1" << endl;
     ofile << "    gray_threshold 0.0" << endl;
     ofile << "    low_error_factor 0.2" << endl;
     ofile << "    minimum_reuse 0.015" << endl;
     ofile << "    nearest_count 10" << endl;
     ofile << "    recursion_limit 20" << endl;
     ofile << "    adc_bailout 0.01" << endl;
     ofile << "    max_sample 0.5" << endl;
     ofile << "    media off" << endl;
     ofile << "    normal off" << endl;
     ofile << "    always_sample 1" << endl;
     ofile << "    pretrace_start 0.08" << endl;
     ofile << "    pretrace_end 0.01" << endl;
     ofile << "}" << endl;
     ofile << "#end" << endl;
     ofile << "}" << endl;
     ofile << endl;

     ofile << "background { color SkyBlue } " << endl << endl;

     JBoundingBox box = mesh->getGeometry()->getBoundingBox();

     Point3D center = box.getCenter();
     double xc = center[0];
     double yc = center[1];
     double zc = center[2];

     ofile << "camera { " << endl;
     ofile << "   location < " << xc << "," << yc << "," << -(box.getLength(1)/tan(M_PI/6) + fabs(zc)+1.0) << ">" <<  endl;
     ofile << "   look_at  < " << xc << "," << yc << "," << -zc << ">" << endl;
     ofile << "}" << endl << endl;

     ofile << "light_source { " << endl;
     ofile << "<" << -xc << "," << yc << "," << -1.5*zc << ">" <<  endl;
     ofile << "  color White " << endl;
     ofile << "} " << endl << endl;

     ofile << "light_source { " << endl;
     ofile << "<" << 0 << "," << yc << "," << -1.5*zc << ">" <<  endl;
     ofile << "color White " << endl;
     ofile << "} " << endl << endl;

     ofile << "light_source { " << endl;
     ofile << "<" << xc << "," << yc << "," << -100.*zc << ">" <<  endl;
     ofile << "color White " << endl;
     ofile << "} " << endl << endl;

     ofile << "light_source { " << endl;
     ofile << "<" << 0 << "," << 0 << "," << -1.5*zc << ">" <<  endl;
     ofile << "color White " << endl;
     ofile << "} " << endl << endl;

     ofile << "light_source { " << endl;
     ofile << "<" << 0.0 << "," << 10*box.getLength(1) << "," << -10*box.getLength(2) << ">" <<  endl;
     ofile << "color White " << endl;
     ofile << "area_light<5,0,0>  <0,0,5> 5, 5 " << endl;
     ofile << "} " << endl << endl;

     ofile << "fog {   " << endl;
     ofile << "     fog_type   2     " << endl;
     ofile << "     distance   50    " << endl;
     ofile << "     color      White " << endl;
     ofile << "     fog_offset 0.1   " << endl;
     ofile << "     fog_alt    2.0   " << endl;
     ofile << "     turbulence 0.8   " << endl;
     ofile << "}" << endl << endl;

     ofile << "plane { " << endl;
     ofile << "  < 0, 1, 0>, " << yc - 0.5*box.getLength(1) << endl;
     ofile << "  pigment { checker color White color Black scale 1.0} " << endl;
     ofile << "  finish  { diffuse 0.9 reflection 0.1 }  " << endl;
     ofile << "} " << endl << endl;

     vector<Point3D> boxedges;
     box.getEdges( boxedges );

     double CylRadius = 0.01*box.getMaxLength();

     ofile << "#declare CylRadius = " << CylRadius << ";" << endl;
     ofile << "#declare MeshBox = object { union { " << endl;

     for( int i = 0; i < 12; i++ ) {
          ofile << "cylinder {<" << boxedges[2*i][0] << "," << boxedges[2*i][1] << "," << -boxedges[2*i][2] << ">,"
          "<" << boxedges[2*i+1][0] << "," << boxedges[2*i+1][1] << "," << -boxedges[2*i+1][2] << ">, CylRadius}" << endl;
     }
     ofile << "}" << endl;
     ofile << " " << "texture { T_Gold_1A } " << endl;
     ofile << "};" << endl;
     ofile <<"MeshBox" << endl;

     double xl = 0.5*box.getMinLength();
     ofile << "#declare XAxis = cylinder {<0,0,0>,<" << xl << ",0, 0>"  << CylRadius << " pigment{color Red} };" << endl;

     ofile << "#declare YAxis = cylinder {<0,0,0>,<0," << xl <<", 0>" << CylRadius << " pigment{color Green}}; "  << endl;

     ofile << "#declare ZAxis = cylinder {<0,0,0>,<0,0," << -xl << ">"  << CylRadius << " pigment{color Blue}}; "   << endl;

     ofile << "#declare Axes = object { union { " << endl;
     ofile << "object{XAxis}" << endl;
     ofile << "object{YAxis}" << endl;
     ofile << "object{ZAxis}" << endl;
     ofile << "}" << endl;
     ofile << "};" << endl;
     ofile << "Axes" << endl;
}


///////////////////////////////////////////////////////////////////////////////

void JMeshPOVExporter :: writemesh( const JMeshPtr &mesh, ofstream &ofile )
{
     ofile <<"#declare SphRadius = 0.001;" << endl;

     ofile << "#declare MeshNodes = object { " << endl;
     ofile << "union {" << endl;
     for( size_t i = 0; i < mesh->getSize(0); i++) {
          const JNodePtr &vertex = mesh->getNodeAt(i);
          const Point3D &xyz = vertex->getXYZCoords();
          ofile << fixed;
          ofile << "sphere {<" << xyz[0] << "," << xyz[1] << "," << -xyz[2] << ">," << "SphRadius}" << endl;
     }
     ofile << "}" << endl;
     ofile << " " << "texture { T_Gold_1A } " << endl;
     ofile << "};" << endl;

     ofile <<"#declare CylRadius = 0.001;" << endl;
     ofile << "#declare MeshEdges = object { " << endl;
///    mesh->getTopology()->collect_edges();
     ofile << "union {  " << endl;
     for( size_t i = 0; i < mesh->getSize(1); i++) {
          JEdgePtr edge = mesh->getEdgeAt(i);
          const Point3D &p0 = edge->getNodeAt(0)->getXYZCoords();
          const Point3D &p1 = edge->getNodeAt(1)->getXYZCoords();
          ofile << "cylinder {<" << p0[0] << "," << p0[1] << "," << -p0[2] << ">,"
          "<" << p1[0] << "," << p1[1] << "," << -p1[2] << ">, CylRadius}" << endl;
     }
     ofile << " } " << endl;
     ofile << "texture{ T_Silver_1A } " << endl;
     ofile << " }; " << endl;

     ofile << "#declare MeshFaces = mesh2 {" << endl;
     ofile << "  vertex_vectors { " << endl;
     ofile << "     " << mesh->getSize(0) << "," << endl;

     for( size_t i = 0; i < mesh->getSize(0); i++) {
          const JNodePtr &vertex = mesh->getNodeAt(i);
          const Point3D &xyz = vertex->getXYZCoords();
          ofile << fixed;
          ofile << "  <" << xyz[0] << "," << xyz[1] << "," << -xyz[2] << ">," << endl;
     }
     ofile << "}" << endl;

     ofile << " face_indices {" << endl;
     ofile << "  " <<  mesh->getSize(2) << "," << endl;

     for( size_t i = 0; i < mesh->getSize(2); i++) {
          JFacePtr face = mesh->getFaceAt(i);
          if( face->isActive() ) {
               ofile << " <";
               ofile << face->getNodeAt(0)->getID() << ",";
               ofile << face->getNodeAt(1)->getID() << ",";
               ofile << face->getNodeAt(2)->getID() << ">," << endl;
          }
     }
     ofile << "  } " << endl;
     ofile << "};" << endl;
}

///////////////////////////////////////////////////////////////////////////////

int
JMeshPOVExporter ::writeFile(const JMeshPtr &mesh, const string &filename)
{
     mesh->pruneAll();

     ofstream ofile( filename.c_str(), ios::out);
     if( ofile.fail() ) {
          cout << "Error: can not open file " << filename << endl;
          return 1;
     }

     writeheader( mesh, ofile );

     ofile << "#include \"model.mesh\" " << endl;

     ofile << "object" << endl;
     ofile << "{" << endl;
     ofile << "    MeshNodes  " << endl;
     ofile << "    texture { T_Stone1 } " << endl;
     ofile << "    rotate<0, 360*clock, 0> " << endl;
     ofile << " } " << endl;
     ofile.close();

     ofile.open("model.mesh", ios::out);
     writemesh( mesh, ofile );

     return 0;
}
