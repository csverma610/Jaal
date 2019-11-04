#include "Mesh.hpp"

#include "BarycentricCoords.hpp"
#include "AllQuadMeshGenerator.hpp"
#include "AllHexMeshGenerator.hpp"


double getVal( double x, double y)
{
    return 0.1*x + 0.3*y;
}


int main()
{

    /*
        int dim[] = {5,5,5};
        double len[] = { 2.0, 2.0, 2.0};
        double org[] = { -1.0, -1.0, -1.0};
        JMeshPtr quadmesh = AllQuadMeshGenerator::getStructuredMesh( dim, len, org );

        JMeshIO::saveAs(quadmesh, "a.xml");
        exit(0);
    */

    JMeshPtr mesh = JMeshIO::readFile("testbary.1.ele");

    vector<Point2D>   cageCoords(4);
    vector<double>    scalarField(4);

    cageCoords[0][0] = 0.0;
    cageCoords[0][1] = 0.0;

    cageCoords[1][0] = 1.0;
    cageCoords[1][1] = 0.0;

    cageCoords[2][0] = 1.00;
    cageCoords[2][1] = 1.00;

    cageCoords[3][0] =  0.0;
    cageCoords[3][1] =  1.0;

    scalarField[0] = getVal(0.0, 0.0);
    scalarField[1] = getVal(1.0, 0.0);
    scalarField[2] = getVal(1.0, 1.0);
    scalarField[3] = getVal(0.0, 1.0);

    JMeanValueCoordinates mvc;

    Point2D qPoint;
    qPoint[0] = 0.5;
    qPoint[1] = 0.5;
    vector<double>  w(4);

    int numnodes = mesh->getSize(0);
    vector<double> shapefunc0, shapefunc2;

    double maxerror = 0.0;
    vector<Vec2D> grad, gradP(4), gradM(4);
    Point2D qPoint2;
    double  h = 0.000001;
    int err;

    for( int i = 0;  i < numnodes; i++) {
        qPoint = mesh->getNodeAt(i)->getXYCoords();
        err = mvc.getCoords(cageCoords, qPoint, shapefunc0);
        cout << shapefunc0[0] << " " << shapefunc0[1] << "  " << shapefunc0[2] << " " << shapefunc0[3] << endl;
        err = mvc.getShapeGradients( cageCoords, qPoint, grad);
/*
        if( !err ) {
            qPoint2[0] = qPoint[0] + h;
            qPoint2[1] = qPoint[1];
            err = mvc.getCoords(cageCoords, qPoint2, shapefunc2);
            maxerror = 0.0;
            if( !err ) {
                cout << "Point2D " << qPoint[0] << " " << qPoint[1] << endl;
                cout << "Analytical " << endl;
                cout << grad[0][0] << "  " << grad[1][0] << " " << grad[2][0] << " " << grad[3][0] << endl;
                cout << grad[0][1] << "  " << grad[1][1] << " " << grad[2][1] << " " << grad[3][1] << endl;
                for( int j = 0; j < 4; j++)  {
                    gradP[j][0] = (shapefunc2[j] - shapefunc0[j] )/h;
                    maxerror = max( maxerror, fabs(grad[j][0]-gradP[j][0]) );
                }
                cout << "Numerical Df/Dx  " << endl;
                cout << gradP[0][0] << "  " << gradP[1][0] << " " << gradP[2][0] << " " << gradP[3][0] << endl;
            }
            qPoint2[0] = qPoint[0];
            qPoint2[1] = qPoint[1]+h;
            err = mvc.getCoords(cageCoords, qPoint2, shapefunc2);
            if( !err) {
                for( int j = 0; j < 4; j++)  {
                    gradP[j][1] = (shapefunc2[j] - shapefunc0[j] )/h;
                    maxerror = max( maxerror, fabs(grad[j][1]-gradP[j][1]) );
                }
                cout << "Numerical Df/Dy  " << endl;
                cout << gradP[0][1] << "  " << gradP[1][1] << " " << gradP[2][1] << " " << gradP[3][1] << endl;
            }
        }
*/
        /*
                double sum = 0.0;
                for( int j = 0; j < 4; j++)
                     sum += shapefunc[j]*scalarField[j];
                double val = getVal(qPoint[0], qPoint[1] );
                double err = fabs(sum - val );
                maxerror = max( maxerror, err);
                cout << qPoint[0] << " " << qPoint[1] << " " << " Calculated " << sum << " Exact " << val << " Error " << err << endl;

        */
    }
    cout << "Maxerror " << maxerror << endl;


    /*
        mvc.setFieldValues( cageCoords, scalarField, mesh, "MVC");
        JMeshVTKExporter mexp;
        mexp.addNodeAttribute("MVC");
        mexp.writeFile(mesh, "mvc.vtk");
    */
}
