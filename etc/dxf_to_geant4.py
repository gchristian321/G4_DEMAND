#!/usr/bin/env python3

import ezdxf
import math
import argparse
import os
from collections import defaultdict

# ---------------------------------------------------------
# Geometry utilities
# ---------------------------------------------------------

def dist(a,b):
    return math.hypot(a[0]-b[0],a[1]-b[1])

def snap(p,tol):
    return (round(p[0]/tol)*tol, round(p[1]/tol)*tol)

def polygon_area(poly):
    s=0
    for i in range(len(poly)):
        x1,y1=poly[i]
        x2,y2=poly[(i+1)%len(poly)]
        s+=x1*y2-x2*y1
    return s/2
    
def same_point(a, b, tol):
    return dist(a, b) <= tol

def is_collinear(a, b, c, tol):
    # 2D cross product magnitude
    cross = (b[0] - a[0])*(c[1] - a[1]) - (b[1] - a[1])*(c[0] - a[0])
    return abs(cross) <= tol

def clean_polygon(poly, tol):
    if len(poly) < 3:
        return poly[:]

    # remove consecutive duplicate points
    out = []
    for p in poly:
        if not out or not same_point(p, out[-1], tol):
            out.append(p)

    # remove duplicate closing point if present
    if len(out) > 1 and same_point(out[0], out[-1], tol):
        out.pop()

    if len(out) < 3:
        return out

    # remove collinear vertices, repeat until stable
    changed = True
    while changed and len(out) >= 3:
        changed = False
        new_out = []
        n = len(out)
        for i in range(n):
            a = out[(i - 1) % n]
            b = out[i]
            c = out[(i + 1) % n]

            # drop b if it is duplicate of neighbor or collinear
            if same_point(a, b, tol) or same_point(b, c, tol) or is_collinear(a, b, c, tol):
                changed = True
                continue

            new_out.append(b)

        out = new_out

    return out

def remove_duplicate_edges(edges,tol):
    seen=set()
    out=[]
    for a,b in edges:
        k=(round(a[0]/tol),round(a[1]/tol),round(b[0]/tol),round(b[1]/tol))
        kr=(k[2],k[3],k[0],k[1])
        if k in seen or kr in seen:
            continue
        seen.add(k)
        out.append((a,b))
    return out

# ---------------------------------------------------------
# Arc discretization
# ---------------------------------------------------------

def arc_points(cx,cy,r,a1,a2,n=48):

    a1=math.radians(a1)
    a2=math.radians(a2)

    while a2<a1:
        a2+=2*math.pi

    pts=[]
    for i in range(n+1):
        t=a1+(a2-a1)*i/n
        pts.append((cx+r*math.cos(t),cy+r*math.sin(t)))
    return pts

# ---------------------------------------------------------
# Read DXF
# ---------------------------------------------------------

def read_dxf_entities(path):

    doc=ezdxf.readfile(path)
    msp=doc.modelspace()

    lines=[]
    circles=[]
    arcs=[]

    for e in msp:

        t=e.dxftype()

        if t=="LINE":
            s=e.dxf.start
            e2=e.dxf.end
            lines.append(((s.x,s.y),(e2.x,e2.y)))

        elif t=="CIRCLE":
            c=e.dxf.center
            circles.append((c.x,c.y,e.dxf.radius))

        elif t=="ARC":
            c=e.dxf.center
            arcs.append((c.x,c.y,e.dxf.radius,
                        e.dxf.start_angle,e.dxf.end_angle))

        elif t=="LWPOLYLINE":

            pts=list(e.get_points())

            for i in range(len(pts)-1):

                x1,y1,_,_,bulge=pts[i]
                x2,y2,_,_,_=pts[i+1]

                if abs(bulge)>1e-12:

                    theta=4*math.atan(bulge)
                    r=dist((x1,y1),(x2,y2))/(2*math.sin(theta/2))

                    segs=arc_points(x1,y1,r,0,theta)

                    for j in range(len(segs)-1):
                        lines.append((segs[j],segs[j+1]))

                else:
                    lines.append(((x1,y1),(x2,y2)))

    return lines,circles,arcs

# ---------------------------------------------------------
# Loop reconstruction
# ---------------------------------------------------------

def build_loops(lines,snap_tol):

    edges=[]
    for a,b in lines:
        a=snap(a,snap_tol)
        b=snap(b,snap_tol)
        if dist(a,b)>1e-9:
            edges.append((a,b))

    edges=remove_duplicate_edges(edges,snap_tol)

    adj=defaultdict(list)

    for i,(a,b) in enumerate(edges):
        adj[a].append((i,b))
        adj[b].append((i,a))

    used=set()
    loops=[]

    for i,(a,b) in enumerate(edges):

        if i in used:
            continue

        used.add(i)

        loop=[a,b]
        prev=a
        cur=b

        while True:

            nxt_edges=[x for x in adj[cur] if x[0] not in used and x[1]!=prev]

            if not nxt_edges:
                break

            idx,nxt=nxt_edges[0]

            used.add(idx)

            loop.append(nxt)

            prev=cur
            cur=nxt

            if cur==loop[0]:
                break

        if len(loop)>3:
            loops.append(loop)

    return loops

# ---------------------------------------------------------
# Write Geant4 C++
# ---------------------------------------------------------

def write_cpp(outfile,outer,holes,circles,thickness):

    basename = os.path.splitext(os.path.basename(outfile))[0]
    func = f"Build{basename}"

    with open(outfile,"w") as f:

        w=f.write

        w('#include "G4ExtrudedSolid.hh"\n')
        w('#include "G4SubtractionSolid.hh"\n')
        w('#include "G4Tubs.hh"\n')
        w('#include "G4LogicalVolume.hh"\n')
        w('#include "G4NistManager.hh"\n\n')

        w(f"G4LogicalVolume* {func}()\n")
        w("{\n")

        w("std::vector<G4TwoVector> outline = {\n")

        for x,y in outer:
            w(f"    {{{x}*mm,{y}*mm}},\n")

        w("};\n\n")

        w(f'G4VSolid* solid = new G4ExtrudedSolid("plate", outline, {thickness}/2*mm, G4TwoVector(),1.0, G4TwoVector(),1.0);\n')

        for i,loop in enumerate(holes):

            w(f"\nstd::vector<G4TwoVector> hole{i} = {{\n")

            for x,y in loop:
                w(f"    {{{x}*mm,{y}*mm}},\n")

            w("};\n")

            w(f'G4VSolid* holeSolid{i} = new G4ExtrudedSolid("hole{i}", hole{i}, {thickness}*mm, G4TwoVector(),1.0, G4TwoVector(),1.0);\n')

            w(f'solid = new G4SubtractionSolid("cut{i}", solid, holeSolid{i});\n')

        for i,(cx,cy,r) in enumerate(circles):

            w(f"""
G4VSolid* circ{i} = new G4Tubs("circ{i}",0,{r}*mm,{thickness}*mm,0,360*deg);
solid = new G4SubtractionSolid("cutc{i}", solid, circ{i}, nullptr, G4ThreeVector({cx}*mm,{cy}*mm,0));
""")

        w("""
auto nist = G4NistManager::Instance();
auto material = nist->FindOrBuildMaterial("G4_Al");

G4LogicalVolume* plateLV =
    new G4LogicalVolume(
        solid,
        material,
        "FrontFrameLV");

return plateLV;
}
""")

# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():

    parser=argparse.ArgumentParser()

    parser.add_argument("dxf")
    parser.add_argument("output_cpp")

    parser.add_argument("--thickness",type=float,required=True)
    parser.add_argument("--hole-max",type=float,default=5)
    parser.add_argument("--snap",type=float,default=1e-4)

    args=parser.parse_args()

    lines,circles,arcs=read_dxf_entities(args.dxf)

    for cx,cy,r,a1,a2 in arcs:
        pts=arc_points(cx,cy,r,a1,a2)
        for i in range(len(pts)-1):
            lines.append((pts[i],pts[i+1]))

    loops=build_loops(lines,args.snap)

    loops.sort(key=lambda L:abs(polygon_area(L)),reverse=True)

    outer=loops[0]
    holes=loops[1:]

    keep_circles=[c for c in circles if c[2]>args.hole_max]

    write_cpp(args.output_cpp,outer,holes,keep_circles,args.thickness)

if __name__=="__main__":
    main()