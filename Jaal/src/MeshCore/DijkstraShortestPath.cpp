#include "StopWatch.hpp"

#ifdef CSV


///////////////////////////////////////////////////////////////////////////////
int JDijkstraShortestPath::atomicOp(LVertex &currnode)
{
    JNodePtr vi = currnode.vertex;
    if (vi == vdst) return 0;

    if( filter ) {
        if( vi != vsrc  ) {
            if( !filter->pass( vi ) ) {
                vdst = vi;
                return 0;
            }
        }
    }

    JNodeSequence vneighs;
    JVertex::getRelations( vi, vneighs );
    int nSize = vneighs.size();
    //   int offset = JMath::random_value((int)0, nSize-1);
    int offset = 0;

    LVertex lv;
    for (int i = 0; i < nSize; i++) {
        JNodePtr vj = vneighs[(i+offset)%nSize];
        if( !vj->isActive() ) continue;
        miter = vmap.find( vj );
        if( miter == vmap.end()) {
            lv.distance = MAXDOUBLE;
            lv.vertex   = vj;
            lv.previous = vi;
            vmap.insert(make_pair(vj,lv));
        } else {
            lv = miter->second;
        }
        assert( lv.vertex == vj );

        double vcost = getCost(currnode, lv);
        if (vcost < lv.distance) {
            lv.previous = vi;
            lv.distance = vcost;
            vmap[vj]    = lv;
            vertexQ.push(lv);
        }
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

void JDijkstraShortestPath::traceback()
{
    nodepath.clear();
    if( vdst == NULL ) return;

    miter = vmap.find(vdst);
    if( miter == vmap.end() ) return;

    LVertex currnode = miter->second;
    while(1) {
        nodepath.push_back( currnode.vertex );
        const JNodePtr &v = currnode.previous;
        if (v == nullptr) break;
        currnode = vmap[v];
    }

    assert( nodepath.front() == vdst );
    assert( nodepath.back()  == vsrc );
}
///////////////////////////////////////////////////////////////////////////////

void JDijkstraShortestPath::fastmarching()
{
    assert( mesh->getAdjTable(0,0) );
    while (!vertexQ.empty()) vertexQ.pop();
    vmap.clear();

    LVertex lv;
    lv.vertex   = vsrc;
    lv.distance = 0.0;
    lv.previous = NULL;
    vmap[vsrc]  = lv;
    vertexQ.push( lv );

    assert( vsrc->isActive() );

    int progress;
    while (!vertexQ.empty()) {
        LVertex currVertex = vertexQ.top();
        vertexQ.pop();
        progress = atomicOp(currVertex);
        if (!progress) break;
    }
}
///////////////////////////////////////////////////////////////////////////////
#endif
