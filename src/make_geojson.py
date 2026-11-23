# 상권 영역 shapefile(EPSG:5181) → 위도·경도 GeoJSON 변환 (geopandas 없이)
# 사용: python make_geojson.py <압축 푼 shapefile 경로(확장자 없이)> <저장할 geojson 경로>
import struct, json, math, sys
a=6378137.0; f=1/298.257222101; e2=f*(2-f); ep2=e2/(1-e2)
lat0=math.radians(38); lon0=math.radians(127); k0=1.0; FE=200000; FN=500000
def M(phi):
    return a*((1-e2/4-3*e2**2/64-5*e2**3/256)*phi-(3*e2/8+3*e2**2/32+45*e2**3/1024)*math.sin(2*phi)+(15*e2**2/256+45*e2**3/1024)*math.sin(4*phi)-(35*e2**3/3072)*math.sin(6*phi))
M0=M(lat0); e1=(1-math.sqrt(1-e2))/(1+math.sqrt(1-e2))
def inv(x,y):
    Mv=M0+(y-FN)/k0; mu=Mv/(a*(1-e2/4-3*e2**2/64-5*e2**3/256))
    p1=mu+(3*e1/2-27*e1**3/32)*math.sin(2*mu)+(21*e1**2/16-55*e1**4/32)*math.sin(4*mu)+(151*e1**3/96)*math.sin(6*mu)+(1097*e1**4/512)*math.sin(8*mu)
    C1=ep2*math.cos(p1)**2; T1=math.tan(p1)**2; N1=a/math.sqrt(1-e2*math.sin(p1)**2); R1=a*(1-e2)/(1-e2*math.sin(p1)**2)**1.5; D=(x-FE)/(N1*k0)
    lat=p1-(N1*math.tan(p1)/R1)*(D**2/2-(5+3*T1+10*C1-4*C1**2-9*ep2)*D**4/24+(61+90*T1+298*C1+45*T1**2-252*ep2-3*C1**2)*D**6/720)
    lon=lon0+(D-(1+2*T1+C1)*D**3/6+(5-2*C1+28*T1-3*C1**2+8*ep2+24*T1**2)*D**5/120)/math.cos(p1)
    return round(math.degrees(lon),6), round(math.degrees(lat),6)
def simplify(pts, tol):
    if len(pts)<5: return pts
    def rdp(s,e,out):
        (x1,y1),(x2,y2)=pts[s],pts[e]; dmax=-1; idx=s
        dx,dy=x2-x1,y2-y1; L=math.hypot(dx,dy)
        for i in range(s+1,e):
            x,y=pts[i]; d=abs(dy*x-dx*y+x2*y1-y2*x1)/L if L else math.hypot(x-x1,y-y1)
            if d>dmax: dmax=d; idx=i
        if dmax>tol: rdp(s,idx,out); rdp(idx,e,out)
        else: out.append(pts[e])
    out=[pts[0]]; rdp(0,len(pts)-1,out); return out if len(out)>=4 else pts
def read(base):
    d=open(base+'.dbf','rb').read(); n,hl,rl=struct.unpack('<IHH',d[4:12]); fields=[]; i=32
    while d[i]!=0x0D:
        fields.append((d[i:i+11].split(b'\0')[0].decode(),d[i+16])); i+=32
    recs=[]
    for r in range(n):
        rec=d[hl+r*rl:hl+(r+1)*rl]; pos=1; o={}
        for name,L in fields: o[name]=rec[pos:pos+L].decode('utf-8','replace').strip(); pos+=L
        recs.append(o)
    s=open(base+'.shp','rb').read(); p=100; geoms=[]
    while p<len(s):
        num,clen=struct.unpack('>ii',s[p:p+8]); c=s[p+8:p+8+clen*2]; p+=8+clen*2
        st=struct.unpack('<i',c[:4])[0]
        if st==0: geoms.append(None); continue
        nparts,npts=struct.unpack('<ii',c[36:44]); parts=list(struct.unpack('<%di'%nparts,c[44:44+4*nparts]))
        xy=struct.unpack('<%dd'%(2*npts),c[44+4*nparts:44+4*nparts+16*npts]); pts=list(zip(xy[0::2],xy[1::2]))
        rings=[pts[parts[k]:(parts[k+1] if k+1<nparts else npts)] for k in range(nparts)]
        geoms.append(rings)
    return recs, geoms
recs, geoms = read(sys.argv[1])
feats=[]
for o,g in zip(recs,geoms):
    if not g: continue
    rings=[[list(inv(x,y)) for x,y in simplify(r,5.0)] for r in g]
    cx,cy=inv(float(o['XCNTS_VALU']),float(o['YDNTS_VALU']))
    feats.append({"type":"Feature","properties":{"상권_코드":int(o['TRDAR_CD']),"상권_코드_명":o['TRDAR_CD_N'],"상권_구분":o['TRDAR_SE_1'],"자치구":o['SIGNGU_CD_'],"행정동":o['ADSTRD_CD_'],"경도":cx,"위도":cy},
                  "geometry":{"type":"Polygon","coordinates":rings}})
json.dump({"type":"FeatureCollection","features":feats},open(sys.argv[2],'w',encoding='utf-8'),ensure_ascii=False,separators=(',',':'))
print(len(feats), inv(197093,453418), inv(200000,500000))
