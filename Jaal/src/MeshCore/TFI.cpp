



void
Jaal::set_tfi_coords(int i, int j, int nx, int ny, JNodeSequence &qnodes)
{
    int offset;

    offset = 0;
    const Point3D &v00 = qnodes[offset]->getXYZCoords();

    offset = i;
    const Point3D &vr0 = qnodes[offset]->getXYZCoords();

    offset = (nx - 1);
    const Point3D &v10 = qnodes[offset]->getXYZCoords();

    offset = j*nx;
    const Point3D &v0s = qnodes[offset]->getXYZCoords();

    offset = j * nx + (nx - 1);
    const Point3D &v1s = qnodes[offset]->getXYZCoords();

    offset = (ny - 1) * nx;
    const Point3D &v01 = qnodes[offset]->getXYZCoords();

    offset = (ny - 1) * nx + i;
    const Point3D &vr1 = qnodes[offset]->getXYZCoords();

    offset = (ny - 1) * nx + (nx - 1);
    const Point3D &v11 = qnodes[offset]->getXYZCoords();

    Point3D vrs;

    double dr = 2.0 / (double) (nx - 1);
    double ds = 2.0 / (double) (ny - 1);

    double r = -1.0 + i*dr;
    double s = -1.0 + j*ds;
    for (int k = 0; k < 3; k++) {
        vrs[k] = TFI::transfinite_blend(r, s,
                                        v00[k], v10[k], v11[k], v01[k],
                                        vr0[k], v1s[k], vr1[k], v0s[k]);
    }
    offset = j * nx + i;
    qnodes[offset]->setXYZCoords(vrs);
}

