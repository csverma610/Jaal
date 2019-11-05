#include "MeshCore/MeshImporter.hpp"
#include "MeshCore/MeshExporter.hpp"
#include <iomanip>

// returns count of non-overlapping occurrences of 'sub' in 'str'
int CountSubString(const std::string& str, const std::string& sub)
{
    if (sub.length() == 0) return 0;
    int count = 0;
    for (size_t offset = str.find(sub); offset != std::string::npos;
            offset = str.find(sub, offset + sub.length()))
    {
        ++count;
    }
    return count;
}

void  JMeshOBJImporter :: readFaceLine(const string &line, vector<size_t> &connect)
{
	/*
    std::vector<std::string> words;
    boost::split(words, line, boost::is_any_of("\t "));
    vector<string> itn;

    connect.clear();
    connect.reserve(words.size());
    for( int i = 0; i < words.size(); i++) {
        if( words[i].empty() )  continue;
        int count = CountSubString(words[i], "/");
        StringTokenizer tok(words[i]);
        itn.clear();
        switch( count)
        {
        case 0:
            connect.push_back(std::stoi( words[i] ));
            break;
        case 1:
            for(StringTokenizer::iterator curTok=tok.begin(); curTok!=tok.end(); ++curTok)
                itn.push_back(*curTok);
            connect.push_back(std::stoi( itn[0]));
            break;
        case 2:
            for(StringTokenizer::iterator curTok=tok.begin(); curTok!=tok.end(); ++curTok)
                itn.push_back(*curTok);
            connect.push_back(std::stoi( itn[0]));
            break;
        }
    }

    int nnodes = connect.size();

//  Since Obj format is based on 1 ...
    for( int i = 0; i < nnodes; i++) 
      connect[i] -= 1;

    switch( nnodes)
    {
    case 3:
        for( int j = 0; j < 3; j++)
            triConnect.push_back(connect[j] );
        break;
    case 4:
        for( int j = 0; j < 4; j++)
            quadConnect.push_back(connect[j] );
        break;
    default:
        if( connect.size() > 4) {
            polyConnect.push_back(connect);
        }
        break;
    }
    */
}

////////////////////////////////////////////////////////////////////////////

JMeshPtr JMeshOBJImporter :: getMesh()
{
    auto mesh = JMesh::newObject();
    size_t index = 0;
    size_t numNodes = nodeCoords.size()/3;

    JNodeSequence connect;

    Point3D p3d;
    JNodeSequence newnodes;
    if( numNodes ) {
        newnodes = JNode::newObjects(numNodes);
        for( size_t i = 0; i < numNodes; i++) {
            p3d[0] = nodeCoords[index++];
            p3d[1] = nodeCoords[index++];
            p3d[2] = nodeCoords[index++];
            newnodes[i]->setXYZCoords(p3d);
            newnodes[i]->setID(i);
        }
        mesh->addObjects(newnodes);
    }

    size_t numTris = triConnect.size()/3;

    index = 0;
    if( numTris ) {
        connect.resize(3);
        auto newtriangles = JTriangle::newObjects(numTris);
        for( size_t i = 0; i < numTris; i++) {
            for( int j = 0; j < 3; j++)
                connect[j] = newnodes[triConnect[index++] ];
            newtriangles[i]->setNodes(connect);
        }
        mesh->addObjects(newtriangles);
    }

    size_t numQuads = quadConnect.size()/4;

    index = 0;
    if( numQuads ) {
        connect.resize(4);
        auto newquads = JQuadrilateral::newObjects(numQuads);
        for( size_t i = 0; i < numQuads; i++) {
            for( int j = 0; j < 4; j++)
                connect[j] = newnodes[quadConnect[index++] ];
            newquads[i]->setNodes(connect);
        }
        mesh->addObjects(newquads);
    }

    size_t numPolys = polyConnect.size();

    index = 0;
    if( numPolys ) {
        auto newpolys = JPolygon::newObjects(numPolys);
        for( size_t i = 0; i < numPolys; i++) {
           int np = polyConnect[i].size();
           connect.resize(np);
            for( int j = 0; j < np; j++)
                connect[j] = newnodes[polyConnect[i][j]];
            newpolys[i]->setNodes(connect);
        }
        mesh->addObjects(newpolys);
    }


    return mesh;
}
////////////////////////////////////////////////////////////////////////////
JMeshPtr JMeshOBJImporter :: getTextureMesh()
{
    JMeshPtr textureMesh = JMesh::newObject();

    size_t index = 0;
    size_t numNodes = texCoords.size()/2;
    Point3D p3d;
    JNodeSequence newnodes, connect;
    if( numNodes ) {
        newnodes = JNode::newObjects(numNodes);
        for( size_t i = 0; i < numNodes; i++) {
            p3d[0] = texCoords[index++];
            p3d[1] = texCoords[index++];
            p3d[2] = 0.0;
            newnodes[i]->setXYZCoords(p3d);
            newnodes[i]->setID(i);
        }
        textureMesh->addObjects(newnodes);
    }

    size_t numTris = triConnect.size()/3;

    index = 0;
    if( numTris ) {
        connect.resize(3);
        auto newtriangles = JTriangle::newObjects(numTris);
        for( size_t i = 0; i < numTris; i++) {
            for( int j = 0; j < 3; j++)
                connect[j] = newnodes[texTriConnect[index++] ];
            newtriangles[i]->setNodes(connect);
        }
        textureMesh->addObjects(newtriangles);
    }

    size_t numQuads = quadConnect.size()/4;

    index = 0;
    if( numQuads ) {
        connect.resize(4);
        auto newquads = JQuadrilateral::newObjects(numQuads);
        for( size_t i = 0; i < numQuads; i++) {
            for( int j = 0; j < 4; j++)
                connect[j] = newnodes[texQuadConnect[index++] ];
            newquads[i]->setNodes(connect);
        }
        textureMesh->addObjects(newquads);
    }
    return textureMesh;

}
////////////////////////////////////////////////////////////////////////////

JMeshPtr JMeshOBJImporter ::readFile(const string &fname)
{
    ifstream infile( fname.c_str(), ios::in);
    if( infile.fail() ) {
        cout << "Warning: Cann't open file " << fname << endl;
        return nullptr;
    }

    //  The codelet is borrowed from TriMesh Software
    JNodePtr vertex;
    JNodeSequence vnodes, connect;
    double  x, y, z;

//    int n1, n2, n3, index = 0;

    nodeCoords.clear();
    texCoords.clear();
    triConnect.clear();
    quadConnect.clear();
    polyConnect.clear();
    texTriConnect.clear();
    texQuadConnect.clear();
    texPolyConnect.clear();

    string str;
    string line;
//    typedef boost::tokenizer<> StringTokenizer;

    size_t numNodes  = 0;
    size_t numFaces  = 0;
    has_texture = 0;

    vector<size_t> faceIndex;
    for(;;) {
        infile >> str;
        if( infile.eof() ) break;
        if( str.compare("#") == 0) getline( infile, line );

        if( str.compare("v") == 0) {
            infile >> x >> y >> z;
            nodeCoords.push_back(x);
            nodeCoords.push_back(y);
            nodeCoords.push_back(z);
            numNodes++;
        }

        if( str.compare("vt") == 0) {
            has_texture = 1;
            infile >> x >> y;
            texCoords.push_back(x);
            texCoords.push_back(y);
        }

        if( str.compare("f") == 0) {
            getline( infile, line );
            readFaceLine( line, faceIndex );
            /*
                        faceIndex.clear();
                        StringTokenizer tok(line);
                        for(StringTokenizer::iterator curTok=tok.begin(); curTok!=tok.end(); ++curTok)
                            faceIndex.push_back(boost::lexical_cast<size_t>(*curTok));
                        readFaceLine( faceIndex );
            */
            numFaces++;
        }
    }

    JMeshPtr mesh = getMesh();
    if( has_texture) {
        JMeshPtr tmesh = getTextureMesh();
        mesh->setAttribute("TextureMesh", tmesh);
    }
    return mesh;
}

//##############################################################################

int
JMeshOBJExporter ::writeFile(const JMeshPtr &mesh, const string &s)
{
    mesh->pruneAll();

    string filename = s;
    ofstream ofile(filename.c_str(), ios::out);
    if( ofile.fail() ) return 1;

    size_t numnodes = mesh->getSize(0);
    size_t numfaces = mesh->getSize(2);
    size_t numCells = mesh->getSize(3);

    for (size_t i = 0; i < numnodes; i++) {
        JNodePtr v = mesh->getNodeAt(i);
        assert( v->isActive() );
        const Point3D &p3d = v->getXYZCoords();
        ofile << fixed << setprecision(precision);
        ofile << "v " << p3d[0] << " " << p3d[1] << " " << p3d[2] << endl;
    }

    JNodeSequence oldConnect, newConnect;
    for (size_t i = 0; i < numfaces; i++) {
        JFacePtr face = mesh->getFaceAt(i);
        if( face->isActive()) {
            if (face->getSize(0) == 4) {
                oldConnect = face->getNodes();
                JQuadrilateral::quad_tessalate(oldConnect, newConnect); // Because of OpenGL
            } else
                newConnect = face->getNodes();

            int nnodes = newConnect.size();
            ofile << "f ";
            for (int j = 0; j < nnodes; j++) {
                size_t vid = newConnect[j]->getID();
                if (vid >= numnodes) {
                    assert(!newConnect[j]->isRemoved());
                    cout << face->getStatus() << endl;
                    cout << "Node indexing out of range " << vid << " Total : " << numnodes << endl;
                    exit(0);
                }
                ofile << vid+1 << " ";
            }
            ofile << endl;
        }
    }

    if( numCells == 0) return 0;

    ofstream volfile;
    string basefile;
    size_t pos1 = filename.rfind(".obj");
    if( pos1 != string::npos) basefile  = filename.substr(0, pos1);
    basefile += "_volmesh.obj";
    volfile.open(basefile.c_str(), ios::out);
    cout  << "An Extended volume mesh is written in :" << basefile << endl;

    for (size_t i = 0; i < numnodes; i++) {
        JNodePtr v = mesh->getNodeAt(i);
        assert( v->isActive() );
        const Point3D &p3d = v->getXYZCoords();
        volfile << "v " << p3d[0] << " " << p3d[1] << " " << p3d[2] << endl;
    }

    for (size_t i = 0; i < numCells; i++) {
        JCellPtr cell = mesh->getCellAt(i);
        int nnodes = cell->getSize(0);
        volfile << "f ";
        for (int j = 0; j < nnodes; j++) {
            volfile << cell->getNodeAt(j)->getID()+1 << " ";
        }
        volfile << endl;
    }

    return 0;
}
/////////////////////////////////////////////////////////////////////////////////
/*
void  JMeshOBJImporter :: readFaceLine( const vector<size_t> &connect)
{
    if( has_texture ) {
        int nnodes = connect.size()/2;
        switch( nnodes)
        {
        case 3:
            for( int j = 0; j < 3; j++) {
                triConnect.push_back(connect[2*j]-1 );
                texTriConnect.push_back(connect[2*j+1]-1 );
            }
            break;
        case 4:
            for( int j = 0; j < 4; j++) {
                quadConnect.push_back(connect[2*j]-1 );
                texQuadConnect.push_back(connect[2*j+1]-1 );
            }
            break;
        default:
            if( connect.size() > 4) {
                polyConnect.push_back( nnodes);
                texPolyConnect.push_back( nnodes);
                for( int j = 0; j < nnodes; j++) {
                    polyConnect.push_back(connect[2*j]-1 );
                    texPolyConnect.push_back(connect[2*j+1]-1 );
                }
            }
            break;
        }
        return;
    }

    int nnodes = connect.size();
    switch( nnodes)
    {
    case 3:
        for( int j = 0; j < 3; j++) {
            triConnect.push_back(connect[j]-1 );
        }
        break;
    case 4:
        for( int j = 0; j < 4; j++) {
            quadConnect.push_back(connect[j]-1 );
        }
        break;
    default:
        if( connect.size() > 4) {
            polyConnect.push_back( nnodes);
            for( int j = 0; j < nnodes; j++) {
                polyConnect.push_back(connect[j]-1 );
            }
        }
        break;
    }
}
*/

