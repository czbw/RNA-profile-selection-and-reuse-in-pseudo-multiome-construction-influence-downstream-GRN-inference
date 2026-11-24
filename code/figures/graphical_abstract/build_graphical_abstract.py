"""Code-rendered AML-GRN graphical abstract. No image-generation model or external art.

Run with the workspace Python runtime. All artwork consists of explicit vector
primitives. Quantitative labels are read from the existing figure manifest.
Network and assignment diagrams are schematic, not inferred network exports.
"""
from pathlib import Path
from math import atan2,cos,sin,pi
import json,html,subprocess,hashlib
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import toColor as HexColor
from PIL import Image

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
MANIFEST=ROOT/'figures/data/manifest.json'
data=json.loads(MANIFEST.read_text(encoding='utf-8'));used={}
def value(key):
    used[key]=data[key]
    return data[key]['value']
V={'regs':value('fig4.set55.n_regulons')}
assert V['regs']==55
W,H=900,360 # 12.5 x 5 inches; exact 5:2 ratio; type is 12-16 pt.
INK='#182B38';MUTED='#53616C';LINE='#D4DCE2';LIGHT='#F2F5F7'
RED='#D11F1F';PURPLE='#5C3696';BLUE='#4292C6';DEEPBLUE='#246991'
FONT=Path('C:/Windows/Fonts')
pdfmetrics.registerFont(TTFont('Arial',str(FONT/'arial.ttf')))
pdfmetrics.registerFont(TTFont('ArialBold',str(FONT/'arialbd.ttf')))
scene=[];text_bounds=[]
def text(x,y,s,size=12,color=INK,bold=False,align='left'):
    face='ArialBold' if bold else 'Arial';width=pdfmetrics.stringWidth(s,face,size)
    left=x-width/2 if align=='center' else x-width if align=='right' else x
    assert left>=8 and left+width<=W-8,(s,left,width)
    assert y-size>=0 and y+3<=H,(s,y)
    scene.append(dict(kind='text',x=x,y=y,s=s,size=size,color=color,bold=bold,align=align))
    text_bounds.append(dict(text=s,left=left,right=left+width,top=y-size,bottom=y+3))
def line(x1,y1,x2,y2,color=INK,width=1.2,dash=False):scene.append(dict(kind='line',x1=x1,y1=y1,x2=x2,y2=y2,color=color,width=width,dash=dash))
def circle(x,y,r,fill='white',stroke=INK,width=1.2):scene.append(dict(kind='circle',x=x,y=y,r=r,fill=fill,stroke=stroke,width=width))
def rect(x,y,w,h,fill='white',stroke=LINE,width=1,radius=0):scene.append(dict(kind='rect',x=x,y=y,w=w,h=h,fill=fill,stroke=stroke,width=width,radius=radius))
def poly(points,fill,stroke=None,width=1):scene.append(dict(kind='poly',points=points,fill=fill,stroke=stroke,width=width))
def arrow(x1,y1,x2,y2,color=MUTED,width=1.25,head=5):
    line(x1,y1,x2,y2,color,width)
    theta=atan2(y2-y1,x2-x1)
    pts=[(x2,y2),(x2-head*cos(theta-pi/6),y2-head*sin(theta-pi/6)),(x2-head*cos(theta+pi/6),y2-head*sin(theta+pi/6))]
    poly(pts,color)
def heading(x,n,s):
    circle(x+9,66,9,INK,INK)
    text(x+9,70,str(n),12,'white',True,'center')
    text(x+26,71,s,14,INK,True)
def network(cx,cy,color,variant=0,scale=1):
    # Explicit generic topology. Node counts/edges have no quantitative meaning.
    pts=[(0,-23),(-27,-7),(27,-7),(-37,22),(0,25),(37,22)]
    edges=[(0,1),(0,2),(1,3),(1,4),(2,5)] if variant==0 else [(0,1),(0,2),(1,4),(2,4),(2,5),(3,4)]
    for a,b in edges:
        x,y=pts[a];xx,yy=pts[b];line(cx+x*scale,cy+y*scale,cx+xx*scale,cy+yy*scale,color,1.35)
    for j,(x,y) in enumerate(pts):circle(cx+x*scale,cy+y*scale,4.6*scale,color if j==0 else 'white',color,1.4)

# Headline and three-stage reading order.
text(450,27,'RNA reuse shapes inferred networks and relapse-associated candidates',16,INK,True,'center')
line(22,40,878,40,LINE,1)
for x in [300,601]:line(x,55,x,306,LINE,1)
heading(22,1,'Pairing changes networks')
heading(321,2,'Test pairing under control')
heading(622,3,'Apply lower-reuse pairing')

# 1: observed pAML discrepancy, using no borrowed biological icons.
text(150,94,'Separately profiled pediatric AML',12,MUTED,False,'center')
text(88,119,'RNA',12,INK,True,'center');text(210,119,'ATAC',12,INK,True,'center')
for x in [68,88,108]:circle(x,136,5.5,'white',MUTED)
for x in [190,210,230]:rect(x-5,131,10,10,'white',MUTED,1.2)
line(88,147,88,155,MUTED);line(210,147,210,155,MUTED);line(88,155,210,155,MUTED)
arrow(149,155,87,174);arrow(149,155,222,174)
text(85,191,'Published link',12,RED,True,'center')
text(222,191,'Metacell',12,PURPLE,True,'center')
network(85,227,RED,0);network(222,227,PURPLE,1)
text(150,273,'Similar activity profiles',13,INK,False,'center')
text(150,292,'Different target sets',14,INK,True,'center')

# 2: two controlled contrasts. Every diagram has RNA left and ATAC right.
# The same eligible RNA pool and ATAC population are drawn in all four cases.
def pairing(cx, y, mapping):
    colors = ['#246991', '#D4833B', '#7A63A5']
    for atac_index, rna_index in enumerate(mapping):
        line(cx-27,y+rna_index*17,cx+27,y+atac_index*17,colors[rna_index],1.5)
    used_rna=set(mapping)
    for j in range(3):
        circle(cx-27,y+j*17,5,colors[j] if j in used_rna else 'white',
               colors[j] if j in used_rna else LINE,1.2)
        rect(cx+22,y+j*17-5,10,10,'white',INK,1.2)
    text(cx-27,y-14,'RNA',12,INK,True,'center')
    text(cx+27,y-14,'ATAC',12,INK,True,'center')

text(451,94,'Same-cell blood and marrow benchmarks',12,MUTED,False,'center')
pairing(386,129,[0,1,2])
pairing(520,129,[1,2,0])
text(451,185,'Partner choice',13,INK,True,'center')
line(336,199,566,199,LINE,0.8)
pairing(386,232,[0,1,2])
pairing(520,232,[1,1,1])
text(451,288,'Reuse ratio',13,INK,True,'center')
text(451,309,'SCENIC+ held fixed',12,MUTED,False,'center')

# 3: schematic network-size contrast, followed by candidate nomination.
# These graphs are not quantitative exports or literal subsets of one another.
def size_network(cx,cy,color,dense=False):
    if dense:
        pts=[(-32,-22),(-9,-31),(18,-27),(39,-12),(-44,4),(-20,-4),
             (4,-9),(28,9),(-34,28),(-8,22),(13,34),(40,29)]
        edges=[(0,1),(0,4),(0,5),(0,6),(1,2),(1,5),(1,6),(1,7),
               (2,3),(2,6),(2,7),(3,7),(3,11),(4,5),(4,8),(4,9),
               (5,6),(5,8),(5,9),(6,7),(6,9),(6,10),(7,10),(7,11),
               (8,9),(9,10),(10,11)]
    else:
        pts=[(0,-22),(-26,-5),(21,-3),(-17,22),(22,21)]
        edges=[(0,1),(0,2),(1,3),(1,4),(2,4)]
    for a,b in edges:
        line(cx+pts[a][0],cy+pts[a][1],cx+pts[b][0],cy+pts[b][1],color,1.05)
    for j,(x,y) in enumerate(pts):
        circle(cx+x,cy+y,4.2,color if j in [0,2] else 'white',color,1.3)

text(750,94,'Longitudinal pediatric AML',12,MUTED,False,'center')
text(682,116,'Published link',12,RED,True,'center')
text(682,133,'Inferred / high reuse',12,MUTED,False,'center')
text(824,116,'Corrected',12,DEEPBLUE,True,'center')
text(824,133,'Random / low reuse',12,MUTED,False,'center')
size_network(682,174,RED,True)
size_network(824,174,DEEPBLUE,False)
# Arrow begins only at the corrected construction.
line(824,205,824,216,DEEPBLUE,1.25)
line(824,216,750,216,DEEPBLUE,1.25)
arrow(750,216,750,226,DEEPBLUE)
text(750,242,f"{V['regs']} relapse-associated",15,DEEPBLUE,True,'center')
text(750,260,'candidate regulons',14,DEEPBLUE,True,'center')
# Adult evidence is stated at module level rather than illustrated as a plot.
text(750,288,'Directionally concordant module',12,INK,False,'center')
text(750,306,'activity changes in adult AML',12,INK,False,'center')

# The principal finding is expressed as a bounded claim, not an accuracy ceiling.
line(22,325,878,325,LINE,1)
text(450,343,'More accurate cell matching does not guarantee better network agreement.',15,INK,True,'center')

# Vector PDF with embedded, selectable Arial text.
c=canvas.Canvas(str(OUT/'graphical_abstract.pdf'),pagesize=(W,H),pageCompression=1,invariant=1)
c.setTitle('AML GRN graphical abstract v3');c.setAuthor('AML-GRN authors; code-assisted draft')
c.setSubject('Original vector schematic rendered from Python; numerical values from the AML-GRN figure manifest')
for o in scene:
    c.saveState();k=o['kind']
    if k=='text':
        c.setFillColor(HexColor(o['color']));c.setFont('ArialBold' if o['bold'] else 'Arial',o['size'])
        fn={'left':c.drawString,'center':c.drawCentredString,'right':c.drawRightString}[o['align']];fn(o['x'],H-o['y'],o['s'])
    elif k=='line':
        c.setStrokeColor(HexColor(o['color']));c.setLineWidth(o['width'])
        if o['dash']:c.setDash(3,3)
        c.line(o['x1'],H-o['y1'],o['x2'],H-o['y2'])
    elif k in ['circle','rect']:
        c.setFillColor(HexColor(o['fill']));c.setStrokeColor(HexColor(o['stroke']));c.setLineWidth(o['width'])
        if k=='circle':c.circle(o['x'],H-o['y'],o['r'],stroke=1,fill=1)
        else:c.roundRect(o['x'],H-o['y']-o['h'],o['w'],o['h'],o['radius'],stroke=1,fill=1)
    else:
        c.setFillColor(HexColor(o['fill']));p=c.beginPath();p.moveTo(o['points'][0][0],H-o['points'][0][1])
        for x,y in o['points'][1:]:p.lineTo(x,H-y)
        p.close();c.drawPath(p,fill=1,stroke=0)
    c.restoreState()
c.showPage();c.save()

esc=html.escape
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="12.5in" height="5in" viewBox="0 0 {W} {H}" role="img" aria-labelledby="title desc">',
'<title id="title">RNA reuse shapes inferred networks and relapse-associated candidates</title>',
'<desc id="desc">Three panels show pairing sensitivity in pediatric AML, controlled same-cell benchmarks, and longitudinal application. RNA is always represented by circles on the left and ATAC by squares on the right. The middle panel compares two partner choices without reuse, then one-to-one pairing with a construction assigning one RNA profile to three ATAC profiles. Unused eligible RNA profiles remain visible as pale circles. The published link network is larger and denser than the corrected network; these are qualitative schematics. An arrow from the corrected network leads to 55 relapse-associated candidate regulons. Text below the candidates states that module activity changes are directionally concordant in adult AML; this does not imply concordance of every individual regulon. Native-pairing references are computational, not experimentally validated networks.</desc>',
f'<rect width="{W}" height="{H}" fill="white"/>']
for o in scene:
    k=o['kind']
    if k=='text':svg.append(f'<text x="{o["x"]}" y="{o["y"]}" fill="{o["color"]}" font-family="Arial, Helvetica, sans-serif" font-size="{o["size"]}" font-weight="{700 if o["bold"] else 400}" text-anchor="{dict(left="start",center="middle",right="end")[o["align"]]}">{esc(o["s"])}</text>')
    elif k=='line':svg.append(f'<line x1="{o["x1"]}" y1="{o["y1"]}" x2="{o["x2"]}" y2="{o["y2"]}" stroke="{o["color"]}" stroke-width="{o["width"]}"/>')
    elif k=='circle':svg.append(f'<circle cx="{o["x"]}" cy="{o["y"]}" r="{o["r"]}" fill="{o["fill"]}" stroke="{o["stroke"]}" stroke-width="{o["width"]}"/>')
    elif k=='rect':svg.append(f'<rect x="{o["x"]}" y="{o["y"]}" width="{o["w"]}" height="{o["h"]}" rx="{o["radius"]}" fill="{o["fill"]}" stroke="{o["stroke"]}" stroke-width="{o["width"]}"/>')
    else:svg.append(f'<polygon points="{" ".join(f"{x},{y}" for x,y in o["points"])}" fill="{o["fill"]}"/>')
svg.append('</svg>');svg='\n'.join(svg);(OUT/'graphical_abstract.svg').write_text(svg,encoding='utf-8')

# Raster formats come from the vector PDF, not an image-generation model.
subprocess.run(['pdftoppm','-singlefile','-r','300','-png',str(OUT/'graphical_abstract.pdf'),str(OUT/'graphical_abstract')],check=True)
im=Image.open(OUT/'graphical_abstract.png').convert('RGB');im.save(OUT/'graphical_abstract.tif',compression='tiff_lzw',dpi=(300,300))
(OUT/'provenance.json').write_text(json.dumps({'created':'2026-11-14','status':'Design draft for author review; not submitted','rendering':'Explicit vector geometry and typography in Python; ReportLab PDF and matching SVG; Poppler raster exports','image_generation_ai_used':False,'llm_assistance':'Codex assisted with design, wording, and rendering code. This is not represented as wholly human-created.','third_party_artwork_used':False,'dimensions_points':[W,H],'aspect_ratio':'5:2','minimum_type_points':12,'manifest':str(MANIFEST),'manifest_sha256':hashlib.sha256(MANIFEST.read_bytes()).hexdigest(),'values':used,'schematic_note':'All graph and assignment glyphs are explanatory illustrations, not representations of measured individual cells or specific inferred networks.','boundaries':['Native pairing is a computational reference, not experimental regulatory ground truth.','Middle-panel diagrams are conceptual controlled contrasts; no quantitative benchmark results are shown.','Network sizes are schematic, not proportional to measured network sizes. A smaller network is not necessarily more accurate.','Adult evidence concerns concordant oriented module activity direction, not concordance of every candidate. Adult RNA results do not validate enhancer edges.']},indent=2),encoding='utf-8')
(OUT/'text_bounds.json').write_text(json.dumps(text_bounds,indent=2),encoding='utf-8')
print('Wrote PDF, SVG, PNG, 300 dpi TIFF and provenance to',OUT)
