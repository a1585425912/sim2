from pathlib import Path
import re,csv,json,copy
from lxml import etree
from docx import Document
from docx.shared import Cm,Pt,RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT,WD_CELL_VERTICAL_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from PIL import Image

BASE=Path(__file__).resolve().parent.parent
WORK=Path(__file__).resolve().parent
OUT=BASE/'实验模型与方法详细报告.docx'
DATA=BASE/'outputs'
M='http://www.w3.org/1998/Math/MathML'
XSL=etree.XSLT(etree.parse(r'C:\Program Files\Microsoft Office\root\Office16\MML2OMML.XSL'))

def el(tag,*children,**attrs):
    e=etree.Element('{'+M+'}'+tag,**attrs)
    for c in children:
        if isinstance(c,str):e.text=(e.text or '')+c
        else:e.append(c)
    return e

GREEK=dict(alpha='α',beta='β',gamma='γ',delta='δ',epsilon='ε',varepsilon='ε',eta='η',theta='θ',lambda_='λ',mu='μ',nu='ν',xi='ξ',rho='ρ',sigma='σ',tau='τ',phi='φ',psi='ψ',omega='ω',Delta='Δ',Sigma='Σ',Lambda='Λ',Omega='Ω',Phi='Φ',chi='χ',pi='π',zeta='ζ')
GREEK['lambda']='λ'
GREEK['varphi']='φ';GREEK['vartheta']='ϑ'
OPS=dict(cdot='·',times='×',otimes='⊗',approx='≈',le='≤',ge='≥',ne='≠',infty='∞',in_='∈',to='→',propto='∝',pm='±',mid='|')

class MathParser:
    def __init__(self,s):self.s=s;self.i=0
    def group(self):
        if self.i<len(self.s) and self.s[self.i]=='{':
            self.i+=1;node=self.parse('}');self.i+=1;return node
        return self.atom()
    def parse(self,stop=''):
        nodes=[]
        while self.i<len(self.s) and (not stop or self.s[self.i]!=stop):
            if self.s[self.i].isspace():self.i+=1;continue
            atom=self.atom();sub=sup=None
            while self.i<len(self.s) and self.s[self.i] in '_^':
                t=self.s[self.i];self.i+=1;v=self.group()
                if t=='_':sub=v
                else:sup=v
            if sub is not None and sup is not None:atom=el('msubsup',atom,sub,sup)
            elif sub is not None:atom=el('msub',atom,sub)
            elif sup is not None:atom=el('msup',atom,sup)
            nodes.append(atom)
        return el('mrow',*nodes)
    def atom(self):
        if self.i>=len(self.s):return el('mrow')
        c=self.s[self.i]
        if c=='{':return self.group()
        if c=='\\':
            self.i+=1
            match=re.match(r'[A-Za-z]+',self.s[self.i:])
            if not match:
                ch=self.s[self.i];self.i+=1
                return el('mo','‖' if ch=='|' else ch)
            cmd=match.group();self.i+=len(cmd)
            if cmd=='frac':return el('mfrac',self.group(),self.group())
            if cmd=='sqrt':return el('msqrt',self.group())
            if cmd in ('hat','dot','ddot','bar'):
                return el('mover',self.group(),el('mo',{'hat':'^','dot':'˙','ddot':'¨','bar':'¯'}[cmd]),accent='true')
            if cmd in ('mathrm','mathbf','mathsf','mathcal','operatorname','text'):
                child=self.group()
                for n in child.iter():
                    if n.tag.endswith('}mi'):n.set('mathvariant','bold' if cmd=='mathbf' else 'normal')
                return child
            if cmd=='begin':
                self.i+=1;end=self.s.index('}',self.i);env=self.s[self.i:end];self.i=end+1
                close='\\end{'+env+'}';end=self.s.index(close,self.i);body=self.s[self.i:end];self.i=end+len(close)
                rows=[el('mtr',*[el('mtd',MathParser(cell).parse()) for cell in row.split('&')]) for row in body.split('\\\\')]
                return el('mrow',el('mo','['),el('mtable',*rows),el('mo',']'))
            if cmd in ('left','right'):return self.atom()
            if cmd in ('sum','int'):return el('mo','∑' if cmd=='sum' else '∫',largeop='true')
            if cmd in ('sin','cos','exp','log','max','min','diag','atan','atan2','Var','Cov','argmin'):
                return el('mi',cmd,mathvariant='normal')
            if cmd in ('quad','qquad'):return el('mspace',width='1em')
            if cmd in GREEK:return el('mi',GREEK[cmd])
            if cmd in OPS:return el('mo',OPS[cmd])
            return el('mi',cmd,mathvariant='normal')
        self.i+=1
        if c.isdigit() or c=='.':
            while self.i<len(self.s) and (self.s[self.i].isdigit() or self.s[self.i]=='.'):c+=self.s[self.i];self.i+=1
            return el('mn',c)
        if c.isalpha():return el('mi',c)
        return el('mo',c)

doc=Document();sec=doc.sections[0]
sec.page_width=Cm(21);sec.page_height=Cm(29.7)
sec.top_margin=Cm(2.1);sec.bottom_margin=Cm(2.0);sec.left_margin=Cm(2.2);sec.right_margin=Cm(2.2)
sec.header_distance=Cm(.9);sec.footer_distance=Cm(.9)

def font(style,name,size,bold=False):
    s=doc.styles[style];s.font.name=name;s.font.size=Pt(size);s.font.bold=bold;s.font.color.rgb=RGBColor(0,0,0)
    rpr=s.element.get_or_add_rPr();f=rpr.find(qn('w:rFonts'))
    if f is None:f=OxmlElement('w:rFonts');rpr.append(f)
    for key in ('ascii','hAnsi','eastAsia'):f.set(qn('w:'+key),name)
    for c in rpr.findall(qn('w:color')):
        for key in ('themeColor','themeShade','themeTint'):c.attrib.pop(qn('w:'+key),None)

for sty,name,size,bold in [('Normal','宋体',11,False),('Title','黑体',22,True),('Subtitle','宋体',11,False),('Heading 1','黑体',15,True),('Heading 2','黑体',12.5,True),('Heading 3','黑体',11,True),('Caption','宋体',9.5,False)]:font(sty,name,size,bold)
normal=doc.styles['Normal'].paragraph_format;normal.line_spacing=1.28;normal.space_after=Pt(6)
normal.widow_control=True
for sty in ('Heading 1','Heading 2','Heading 3'):
    pf=doc.styles[sty].paragraph_format;pf.space_before=Pt(13);pf.space_after=Pt(7);pf.keep_with_next=True
doc.styles['Caption'].paragraph_format.space_after=Pt(8)
hp=sec.header.paragraphs[0];hp.text='动态浮标联合跟踪UAV与相干合成实验说明';hp.alignment=WD_ALIGN_PARAGRAPH.RIGHT
for r in hp.runs:r.font.size=Pt(8);r.font.color.rgb=RGBColor(0,0,0)
fp=sec.footer.paragraphs[0];fp.alignment=WD_ALIGN_PARAGRAPH.CENTER
fp.add_run('第 ')
field=OxmlElement('w:fldSimple');field.set(qn('w:instr'),'PAGE');fp._p.append(field)
fp.add_run(' 页')
for r in fp.runs:r.font.size=Pt(9)

eqcount=0
def p(text='',bold=False):
    par=doc.add_paragraph();run=par.add_run(text);run.bold=bold;return par
def h(text,level=1):return doc.add_heading(text,level)
def eq(s,explain=None):
    global eqcount
    eqcount+=1
    math=el('math',MathParser(s).parse(),display='block')
    omml=XSL(math).getroot()
    par=doc.add_paragraph();par.alignment=WD_ALIGN_PARAGRAPH.CENTER
    par.paragraph_format.space_before=Pt(4);par.paragraph_format.space_after=Pt(7)
    par._p.append(copy.deepcopy(omml));r=par.add_run('   （'+str(eqcount)+'）');r.font.size=Pt(9)
    if explain:p(explain)

def table(headers,rows,widths=None):
    t=doc.add_table(rows=1,cols=len(headers));t.alignment=WD_TABLE_ALIGNMENT.CENTER;t.autofit=False
    if widths is None:widths=[16.6/len(headers)]*len(headers)
    for i,w in enumerate(widths):t.columns[i].width=Cm(w)
    for i,x in enumerate(headers):t.rows[0].cells[i].text=str(x)
    for row in rows:
        cells=t.add_row().cells
        for i,x in enumerate(row):cells[i].text=str(x)
    for ri,row in enumerate(t.rows):
        prop=row._tr.get_or_add_trPr();nsplit=OxmlElement('w:cantSplit');prop.append(nsplit)
        if ri==0:repeat=OxmlElement('w:tblHeader');prop.append(repeat)
        for ci,cell in enumerate(row.cells):
            cell.width=Cm(widths[ci]);cell.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
            pr=cell._tc.get_or_add_tcPr();borders=OxmlElement('w:tcBorders')
            for side in ('top','left','bottom','right'):
                item=OxmlElement('w:'+side);item.set(qn('w:val'),'single');item.set(qn('w:sz'),'5');item.set(qn('w:color'),'D9D9D9');borders.append(item)
            pr.append(borders);marg=OxmlElement('w:tcMar')
            for side,n in [('top',85),('bottom',85),('left',95),('right',95)]:
                item=OxmlElement('w:'+side);item.set(qn('w:w'),str(n));item.set(qn('w:type'),'dxa');marg.append(item)
            pr.append(marg)
            if ri==0 or ri%2==0:
                shade=OxmlElement('w:shd');shade.set(qn('w:fill'),'DCE6F1' if ri==0 else 'F5F7FA');pr.append(shade)
            for par in cell.paragraphs:
                par.paragraph_format.space_after=Pt(1);par.paragraph_format.line_spacing=1.12
                par.alignment=WD_ALIGN_PARAGRAPH.LEFT if ci==0 or len(str(cell.text))>25 else WD_ALIGN_PARAGRAPH.CENTER
                for run in par.runs:run.font.size=Pt(9);run.bold=ri==0
    doc.add_paragraph().paragraph_format.space_after=Pt(1)
    return t

figcount=0
def pic(path,caption,width=16.2):
    global figcount
    figcount+=1
    par=doc.add_paragraph();par.alignment=WD_ALIGN_PARAGRAPH.CENTER;par.paragraph_format.keep_with_next=True
    run=par.add_run();run.add_picture(str(path),width=Cm(width))
    for node in run._r.xpath('.//wp:docPr'):node.set('descr',caption)
    cap=doc.add_paragraph('图'+str(figcount)+'  '+caption,style='Caption');cap.alignment=WD_ALIGN_PARAGRAPH.CENTER

summary=list(csv.DictReader((DATA/'summary.csv').open(encoding='utf-8-sig')))
def get(kind,method='GP_Doppler_q1',hs=1.5,timing='measured_cpu',freq=300,phase='all'):
    return next(x for x in summary if float(x['hs'])==hs and x['trajectory']==kind and x['method']==method and x['timing']==timing and float(x['frequency_mhz'])==freq and x['phase']==phase)

doc.add_paragraph('FMCW距离率辅助动态浮标联合跟踪UAV实验',style='Title')
doc.add_paragraph('模型 方法 数学推导 完整实验结果与实现边界',style='Subtitle')
p('版本日期 2026年11月1日    实验状态 完整批次已结束')
p('本报告面向需要理解并继续研究该仿真的项目使用者，说明这次实验到底模拟了什么、估计了什么，以及新增FMCW径向速度如何进入轨迹估计和相干发射。内容依据当前实际执行源码及最终20种子统计编写，不把尚未实现的设想写成已完成算法。')
p('主要结论：距离率紧耦合在悬停、圆周和混合机动中明显改善速度、波动、未来几何与相干收益；但长直线和爬升的远段发生测速异常与轨迹跳变，固定高斯观测模型的置信度失配。当前版本尚未同时满足全程可靠、200毫秒帧内时序与1 GHz高相干要求。',True)
table(['项目','本轮实际内容'],[
('完整规模','20种子 × 2海况 × 5轨迹 × 5方法 全部0至90秒'),
('结果文件','32000行跟踪与相干指标 1600行测速指标'),
('物理仿真','连续双时间尺度浮标运动 原始导航传感器 新增复数FMCW回波'),
('主估计器','持久联合高斯到达因子 WNOJ先验 三次迭代MAP更新'),
('本轮新增','相对径向速度因子及移动浮标速度联合约束'),
('未采用','本轮没有重跑IMM或CI 没有Student t观测模型或功率优化')],[3,13.6])
h('阅读顺序',2)
p('第一至三章定义任务、坐标和运动真值；第四至六章介绍原始观测、导航与FMCW测速；第七至九章推导速度因子、GP动力学和联合求解；第十至十二章说明预测、相干评价、完整实验和统计；最后给出结果、限制、源码映射及正式参考文献。全部公式均为可编辑Word数学对象。')

h('一 研究任务与模型层次')
p('最终任务有两个相互关联而不等价的目标：一是获得误差较小、波动受控且不确定度可信的UAV轨迹；二是让各浮标利用估计后的未来相对几何设置发射相位，在UAV位置形成较高接收功率。位置RMSE较小不自动保证相干增益最大，平均接收功率增加也不自动证明全程轨迹可靠。')
table(['层次','输入与输出','是否访问运动真值'],[
('前向物理与传感器','海场和轨道生成原始传感器与ADC','可以 只负责生成与评价'),
('导航参考','伪距 IMU 地磁生成误差坐标参考','不可以'),
('联合估计','原始观测与到达先验生成目标及浮标后验','不可以'),
('未来预测与命令','后验生成未来双方位置和传播距离命令','不可以'),
('仿真评价','真实未来距离生成实际复场并比较误差','可以 结果不反馈估计')],[3,8.6,5])
p('“不是给真值加噪声的假导航”具体指：浮标导航位置和速度不是由真值直接扰动产生，而是从伪距、比力、角速度、地磁读数解算。新增距离率也来自复数回波的频谱处理。原角距观测仍采用几何前向函数加原始测量噪声；这属于传感器级抽象，不是已经完成全部雷达射频与天线阵列仿真。')
p('本轮是集中联合信息参考。20个浮标的数据进入一个联合估计器，再用于各节点预测；尚未实现与有限轮CI具有同等通信约束的分布式GP。因此125毫秒通信只是沿用的时延假设，不是已完成网络部署的证据。')

h('二 坐标 状态和符号')
p('采用局部右手XYZ世界坐标，Z轴向上，长度单位为米，时间单位为秒。R_i把浮标机体系向量旋转到世界系，重力取g=[0,0,−9.81]ᵀ。所有导航传感器、雷达参考点与等效通信发射点共位，杆臂为零；浮标姿态仍随海浪变化。')
table(['符号','含义'],[
('p_u v_u a_u','UAV世界系位置 速度 加速度 各为三维'),
('p_b,i v_b,i R_i','第i浮标位置 速度 机体到世界旋转'),
('c_i s_i','浮标慢漂移中心与波频平移 总位置为二者之和'),
('r_i u_i','浮标至UAV距离及该方向的单位视线'),
('q 或 q_j','jerk功率谱密度 单位m²/s⁵ 不是速度'),
('λ_r λ_c','FMCW雷达波长与通信载波波长 二者分开'),
('P H R_z','状态协方差 观测雅可比 测量噪声协方差'),
('χ U N','联合未知量 目标克隆维数 浮标数20')],[4.2,12.4])
eq(r'\mathbf{x}_u=[\mathbf{p}_u^T,\mathbf{v}_u^T,\mathbf{a}_u^T]^T')
eq(r'\delta\mathbf{n}_i=[\delta\mathbf{p}_i^T,\delta\mathbf{v}_i^T,\delta\theta_i^T,\delta\mathbf{b}_{a,i}^T,\delta\mathbf{b}_{g,i}^T,\delta c_i]^T')
p('单浮标导航误差为16维，最后一维是已乘光速的GNSS接收机钟差，单位米。右姿态误差定义为R_true=R_ref Exp([δθ]×)，其中[·]×是叉乘矩阵。钟差符号c与慢中心c_i在不同上下文中使用，源码中clock明确区分。')

h('三 连续运动真值模型')
h('3 1 UAV的五类运动',2)
p('UAV真值仅用于生成观测和评价。估计器不会获得“当前悬停”“当前加速”等真实标签，也不会获得轨迹切换时刻。所有方法使用相同轨迹。')
table(['轨迹','实际定义'],[
('悬停','固定位置[0 0 20]米 世界速度为零'),
('直线','p=[15t −200 120]米 速度[15 0 0]米每秒'),
('圆周','p=[200sin(0.075t) −200cos(0.075t) 120]米 速度大小15米每秒'),
('混合','初始[0 −200 120]米 15至20秒加速 35至55秒转弯 65至70秒刹停'),
('爬升','直线水平运动基础上 20至22秒建立2米每秒垂直速度 40至42秒停止爬升')],[2.4,14.2])
eq(r's(t)=3\min(\max(t-15,0),5)-3\min(\max(t-65,0),5)')
eq(r'\psi(t)=0.1\min(\max(t-35,0),20),\quad \mathbf{v}(t)=[s\cos\psi,s\sin\psi,0]^T')
p('混合位置由上述速度在1毫秒网格上累计梯形积分，再对查询时刻插值。此网格只定义前向真值，不是估计器的轨迹拟合输入。爬升的垂直速度为clip(t−20,0,2)−clip(t−40,0,2)，相应位置通过分段二次积分构造。')
h('3 2 浮标慢漂移',2)
eq(r'\mathbf{p}_{b,i}(t)=\mathbf{c}_i(t)+\mathbf{s}_i(t)')
eq(r'\mathbf{v}_{c,i}(t)=\bar{\mathbf{v}}+\sum_{l=1}^{3}[\mathbf{A}_{il}\cos(\nu_l t)+\mathbf{B}_{il}\sin(\nu_l t)]')
eq(r'\mathbf{c}_i(t)=\mathbf{c}_i(0)+\bar{\mathbf{v}}t+\sum_{l=1}^{3}[\frac{\mathbf{A}_{il}}{\nu_l}\sin(\nu_l t)+\frac{\mathbf{B}_{il}}{\nu_l}(1-\cos(\nu_l t))]')
p('第二式对时间积分即得第三式，积分常数保证t=0时回到部署中心。平均流速为[0.3,0.1,0]米每秒，三个周期为180、120、80秒。A与B只在每个种子开始时抽取一次，水平分量标准差为0.05/√3米每秒，垂直慢速度没有额外随机谱项。因此这是平滑有限谱随机漂移，不是每个时刻独立跳动，也不是OU过程。')
p('20个慢中心初始部署在水平半径100米的圆周上，高度1.5米。相同种子在两个海况中保留一致随机结构，海况主要改变波浪幅度。')
h('3 3 共享方向谱与波频运动',2)
eq(r'\eta(\mathbf{r},t)=\sum_{m=1}^{48}A_m\cos(\omega_m t-\mathbf{k}_m^T\mathbf{r}+\phi_m)')
eq(r'\sum_m\frac{A_m^2}{2}=\frac{H_s^2}{16},\qquad \omega_m^2=gk_m\tanh(k_m d)')
p('谱形采用JONSWAP形式，48个分量，峰值周期6秒，谱峰因子3.3，水深50米，波向中心30度、角度离散标准差30度。频率范围为峰值角频率的0.6至3倍。色散方程通过Newton迭代求波数。Hs为波面有效波高，不能直接说浮标升沉标准差必然等于Hs/4，因为浮标还经过响应函数。海洋低频与波频区分的研究背景见文献[1]。')
eq(r'\mathbf{s}_i(t)=\operatorname{Re}\sum_m\mathbf{H}_{p,m}\exp(j\psi_{im}),\quad \psi_{im}=\omega_m t-\mathbf{k}_m^T\mathbf{c}_{i,xy}(t)+\phi_m')
eq(r'H_z(\omega)=\frac{1}{1-(\omega/\omega_n)^2+2j\zeta\omega/\omega_n}')
p('升沉自然周期8秒，阻尼比0.2；水平响应取0.7倍谱幅并采用正交相位。具体响应为Hp,m=[−j0.7Am cosβm,−j0.7Am sinβm,AmHz(ωm)]。这些是明确的轻量响应假设，没有实测浮标水动力标定，不能借文献船型参数声称已获得真实浮标响应。')
eq(r'\dot\psi_{im}=\omega_m-\mathbf{k}_m^T\mathbf{v}_{c,i,xy},\qquad \ddot\psi_{im}=-\mathbf{k}_m^T\mathbf{a}_{c,i,xy}')
eq(r'\dot{\mathbf{s}}_i=\operatorname{Re}\sum_m j\dot\psi_{im}\mathbf{H}_{p,m}e^{j\psi_{im}}')
eq(r'\ddot{\mathbf{s}}_i=\operatorname{Re}\sum_m[-\dot\psi_{im}^{2}+j\ddot\psi_{im}]\mathbf{H}_{p,m}e^{j\psi_{im}}')
p('这三式由链式法则连续求导而来，使位置、速度、加速度属于同一条连续轨迹；同一波场在不同浮标慢中心处采样，保留空间相关性。没有把真值慢快分量提供给后续预测器。')
h('3 4 浮标姿态与零杆臂',2)
p('姿态角也由共享波场的简化复响应产生，滚转与俯仰与波面斜率相关，偏航含随机初始航向及较小波频项。欧拉角构造R=Rz(yaw)Ry(pitch)Rx(roll)，再由角度导数计算机体系角速度。')
eq(r'\omega_B=[\dot\varphi-\dot\psi\sin\vartheta,\dot\vartheta\cos\varphi+\dot\psi\sin\varphi\cos\vartheta,-\dot\vartheta\sin\varphi+\dot\psi\cos\varphi\cos\vartheta]^T')
p('这里φ、ϑ、ψ分别为滚转、俯仰和偏航。关系保证Ṙ=R[ωB]×。零杆臂意味着绕参考点的旋转不产生额外ω×lever平移速度；不意味着姿态恒定，也不意味着姿态误差可以从IMU、磁场或雷达角度模型中删掉。')

h('四 原始传感器与导航参考')
h('4 1 原始测量方程',2)
eq(r'z_{\rho,ij}=\|\mathbf{s}_j-\mathbf{p}_{b,i}\|+c_{clock,i}+n_{\rho,ij}')
eq(r'\mathbf{f}_m=R_i^T(\mathbf{a}_{b,i}-\mathbf{g})+\mathbf{b}_{a,i}+\mathbf{n}_a,\quad \omega_m=\omega_B+\mathbf{b}_{g,i}+\mathbf{n}_g')
eq(r'\mathbf{m}_m=R_i^T\mathbf{m}_{world}+\mathbf{n}_m')
p('GNSS为八颗固定合成可见卫星的伪距，而非已知真实浮标坐标。钟差真值按离散OU过程生成，估计模型用随机游走近似；未模拟大气、多径和载波模糊度。IMU以50 Hz输出区间末可用的比力和角速度，用中点真值近似区间平均。地磁为已知非竖直向量[20,0,−45]微特斯拉，不输入真航向。')
table(['参数','设置'],[
('伪距与角距更新','5 Hz'),('伪距标准差','0.30 m'),('加速度白噪声密度','0.015 m/s²/√Hz'),('陀螺白噪声密度','0.0003 rad/s/√Hz'),('加速度偏置初始标准差','0.03 m/s²'),('陀螺偏置初始标准差','0.001 rad/s'),('两类偏置随机游走密度','0.0003 m/s²/√s 与 0.00001 rad/s/√s'),('地磁单分量标准差','0.15 μT'),('原雷达距离与角度标准差','0.8 m 与 0.35°')],[7,9.6])
p('偏置随机游走按b(k+1)=b(k)+σrw√Δt ξ生成。连续白噪声密度为σ时，区间平均读数的标准差是σ/√Δt；不能把噪声密度直接当作每个采样点的标准差。')
eq(r'\mathbf{q}_B=R_i^T(\mathbf{p}_u-\mathbf{p}_{b,i}),\quad h_{radar}=[\|\mathbf{q}_B\|,\atan2(q_y,q_x),\atan2(q_z,\sqrt{q_x^2+q_y^2})]^T')
p('原角距通道按上述原始球坐标函数生成，并叠加原测量噪声。角度残差归一到[−π,π)。新FMCW通道只补充速度；不能把这轮结果称为已经从同一ADC同时解出了全部角度和亚米级距离。')
h('4 2 从观测初始化',2)
eq(r'(\hat{\mathbf{p}}_0,\hat c_0)=\argmin_{\mathbf{p},c}\sum_j\frac{[z_{\rho,j}-\|\mathbf{s}_j-\mathbf{p}\|-c]^2}{\sigma_\rho^2}')
p('位置与钟差从局部原点和零钟差开始，做六次Gauss–Newton迭代。初始倾角与航向通过首个加速度及磁向量构造正交基得到，并给予10度姿态标准差以容纳波浪加速度的影响。初速度和两类偏置均初始化为零，速度先验单轴标准差1 m/s。没有使用真实位置、速度、姿态或偏置。首帧已用于初始化的GNSS与地磁不会再次加入联合更新。')
h('4 3 误差状态惯性递推',2)
eq(r'\dot{\mathbf{p}}=\mathbf{v},\quad \dot{\mathbf{v}}=R(\mathbf{f}_m-\mathbf{b}_a)+\mathbf{g},\quad \dot R=R[\omega_m-\mathbf{b}_g]_\times')
eq(r'\delta\dot{\mathbf{v}}=-R[\mathbf{f}]_\times\delta\theta-R\delta\mathbf{b}_a-R\mathbf{n}_a')
eq(r'\delta\dot\theta=-[\omega]_\times\delta\theta-\delta\mathbf{b}_g-\mathbf{n}_g,\qquad F_n\approx I+A\Delta t+\frac{A^2\Delta t^2}{2}')
p('第二、三式由右姿态小扰动和一阶线性化得到。平移用中点姿态推进；过程噪声保留加速度对白噪声位置与速度的耦合：Qpp=σa²Δt³/3、Qpv=σa²Δt²/2、Qvv=σa²Δt，姿态为σg²Δt，偏置与钟差为随机游走密度平方乘Δt。这是明确的有限阶离散近似，不是全16维精确离散化。')
p('参考ESKF使用Joseph协方差更新，并在姿态注入后进行一阶误差重置。其作用是给联合估计提供坐标参考和原始惯性转移；当前ESKF后验不会作为一份新的独立导航观测重复输入。文献[5]给出基本线性高斯更新背景，但本项目非线性小角近似仍需经验检验。')

h('五 FMCW回波模型和测速方法')
h('5 1 完整扫频与采样窗口',2)
p('77 GHz是ADC起点的瞬时频率。每个发射线性段持续50 μs，前10 μs为不记录的前置段，末40 μs为ADC采样窗。5 MHz指采样窗内有效扫频带宽，因此斜率S=5 MHz/40 μs，完整TX带宽为6.25 MHz。源码历史字段chirp_s和部分注释写40 μs，其数值实际用于ADC窗长度；不能把它误解成没有前置段的完整TX时长。')
table(['量','数值与含义'],[
('名义频率与波长','77 GHz  λr约3.8934 mm'),('完整扫频与重复周期','50 μs  无空隙的理想锯齿线性段'),('采样窗与前置段','40 μs 与10 μs'),('单窗快时间采样','64点  fs=1.6 MHz'),('有效及完整发射带宽','5 MHz 与6.25 MHz'),('每CPI扫频数','256'),('名义CPI与半支持区间','12.8 ms 与6.395 ms'),('慢时间FFT','512点 即2倍零填充'),('参考单采样SNR','250 m处5 dB 幅度随距离平方衰减'),('融合中的固定σD','0.15 m/s 不是FFT理论精度')],[7,9.6])
h('5 2 从延迟信号推导解调相位',2)
eq(r's_m(\tau)=\exp\{j\phi_m+j2\pi[f_c\tau+\frac{S\tau^2}{2}]\},\quad \delta=\frac{2r}{c}')
eq(r's_m(\tau)s_m^{*}(\tau-\delta)=\exp\{j2\pi[(f_c+S\tau)\delta-\frac{S\delta^2}{2}]\}')
p('第二式由两项相位相减得到，段初相φm抵消。雷达回波是往返传播，因此延迟含2r/c。对于40 μs采样窗内的τ，只要τ−δ仍在同一发射线性段，式子成立；固定10 μs前置段保证本次所有目标距离均满足同段接收条件。跨扫频起始相位可选为连续，但同段混频结果不依赖这些常数。')
eq(r'r_{imn}=r_{im}+\dot r_{im}t_n,\quad \delta_{imn}=\frac{2r_{imn}}{c}')
eq(r'y_{imn}=A_0(\frac{r_0}{r_{imn}})^2\exp\{j2\pi[(f_c+S\tau_n)\delta_{imn}-\frac{S\delta_{imn}^2}{2}]\}+w_{imn}')
p('m为扫频，n为快时间采样。每个扫频中心的距离及距离率由同一连续双方轨迹生成；仅在40 μs窗内一阶展开距离。扫频之间保留真实机动、浮标升沉和加速度造成的相位演化。热噪声加在复数I/Q采样中，不在距离率真值上直接加一个人为速度噪声。该模型忽略光行时内运动的高阶项和单窗加速度二阶项，不是精确电磁或射频电路模拟。')
eq(r'\operatorname{SNR}(r)=5+40\log_{10}(\frac{250}{r})\quad \mathrm{dB}')
p('幅度按r⁻²变化，所以功率按r⁻⁴变化，转换为dB即产生40log10项。复噪声功率为10^(−5/10)，I、Q各占一半。参考SNR为批次前固定设计假设，未通过实际天线增益、UAV RCS、噪声系数和发射功率闭合链路预算，不能据此宣称某台77 GHz硬件能在1.5 km可靠探测；本实验1 W限制属于通信而非新增感知雷达。')
h('5 3 距离多普勒处理与全图选峰',2)
p('流程为快慢两轴Hann加窗、快时间64点FFT、慢时间512点FFT并居中，然后对每个浮标的整个二维功率图取最大峰。既不用真实距离设置搜索窗，也不用真实速度挑峰；不做CFAR、低质量帧剔除或测量门控。它是已知单目标场景中的基础选峰器，不代表已经解决复杂海面目标检测与关联。')
eq(r'\epsilon=\frac{L_{-1}-L_{+1}}{2(L_{-1}-2L_0+L_{+1})},\quad \hat f_D=\frac{k_D-N_D/2+\epsilon}{N_DT_r}')
eq(r'\hat{\dot r}=\frac{\lambda_r\hat f_D}{2},\quad \Delta v=\frac{\lambda_r}{2MT_r},\quad |v|_{max}=\frac{\lambda_r}{4T_r}')
p('L为峰所在距离单元内相邻三个慢时间对数功率。抛物线顶点给出亚栅格偏移ε。物理速度FFT间隔约0.1521 m/s；零填充后的网格约0.0760 m/s，并不把物理分辨率翻倍。无模糊速度约±19.467 m/s；高SNR下估计精度可小于栅格间距，但不能把分辨率或零填充网格直接当作测量标准差。')
p('快时间复I/Q采用已知正拍频的[0,fs)映射，不做快轴fftshift，距离周期约1918.672 m；它不同于通常有符号频带下的正距离上界，也不同于前置段同扫频接收范围1498.962 m。全部200条轨迹的ADC期间距离上界为1450.552 m，径向速度上界为15.583 m/s，均在本次模型域内；核对只用于物理有效性，不参与选峰或删帧。')
pic(DATA/'range_doppler.png','复数ADC产生的完整距离多普勒功率图 示例为首种子浮标1 t等于20秒')

h('六 相对径向速度与几何信息')
h('6 1 从距离求导到速度因子',2)
eq(r'\mathbf{q}_i=\mathbf{p}_u-\mathbf{p}_{b,i},\quad r_i=\sqrt{\mathbf{q}_i^T\mathbf{q}_i},\quad \mathbf{u}_i=\frac{\mathbf{q}_i}{r_i}')
eq(r'\dot r_i=\frac{\mathbf{q}_i^T\dot{\mathbf{q}}_i}{r_i}=\mathbf{u}_i^T(\mathbf{v}_u-\mathbf{v}_{b,i})')
p('这是核心物理关系。测到的是UAV相对浮标的速度投影，不是UAV世界速度。若UAV悬停，仍有ṙ=−uᵀvb，浮标慢漂移和海浪升沉会进入测速。零杆臂时没有额外转动平移项，但姿态误差仍通过导航和角距观测耦合。')
eq(r'e_{D,i}=z_{D,i}-\mathbf{u}_i^T(\mathbf{v}_u-\mathbf{v}_{b,i}),\quad J_D=\frac12\sum_i\frac{e_{D,i}^2}{\sigma_{D,i}^2}')
p('它首先是新增观测似然，不是运动先验。固定σD=0.15 m/s为探索性名义尺度，未用真值校准。本轮角距与新FMCW热噪声按独立通道生成，因此可以这样叠加；若未来改成同一回波同时提取角距与速度，就应重新考虑联合观测协方差和相关性。')
h('6 2 观测雅可比逐步推导',2)
eq(r'd\mathbf{u}_i=\frac{I-\mathbf{u}_i\mathbf{u}_i^T}{r_i}d\mathbf{q}_i,\quad \mathbf{g}_i=\frac{(I-\mathbf{u}_i\mathbf{u}_i^T)(\mathbf{v}_u-\mathbf{v}_{b,i})}{r_i}')
eq(r'd\dot r_i=\mathbf{g}_i^T(d\mathbf{p}_u-d\mathbf{p}_{b,i})+\mathbf{u}_i^T(d\mathbf{v}_u-d\mathbf{v}_{b,i})')
eq(r'H_{D,i}=[\mathbf{g}_i^T,\mathbf{u}_i^T,0,-\mathbf{g}_i^T,-\mathbf{u}_i^T]')
p('先对单位向量q/||q||应用商法则，得到投影矩阵I−uuᵀ；再对uᵀ(vu−vb)用乘积法则。H按[pu,vu,au,pb,vb]排列。负号说明同一径向测量同时约束两方，不能先扣一个浮标速度估计值，再把它当成完全已知。姿态和偏置通过联合协方差间接被更新。')
h('6 3 多浮标为何能补充三维速度',2)
eq(r'\mathbf{y}=U\mathbf{v}_u+\mathbf{n},\quad U=[\mathbf{u}_1,\ldots,\mathbf{u}_N]^T')
eq(r'\hat{\mathbf{v}}_u=(U^TR_z^{-1}U)^{-1}U^TR_z^{-1}\mathbf{y},\quad P_v=(U^TR_z^{-1}U)^{-1}')
p('这里为了说明几何，暂假设视线和浮标速度已知，并令yi=zD,i+uiᵀvb,i。单站只有一个速度投影，多站方向互补时U可满列秩；远距离视线趋于平行时，信息矩阵的某些特征值变小，横向或垂直速度更敏感。真实联合模型还包含导航不确定度，因此不能直接把上述条件协方差当最终置信度。文献[2]支持距离与距离率联合提供运动信息，但其双基地五参数场景不直接证明本项目任意三轴运动都可观。')

h('七 UAV轨迹先验和联合状态')
h('7 1 白噪声jerk模型',2)
eq(r'\dot{\mathbf{p}}=\mathbf{v},\quad \dot{\mathbf{v}}=\mathbf{a},\quad \dot{\mathbf{a}}=\mathbf{w},\quad E[\mathbf{w}(t)\mathbf{w}(s)^T]=qI\delta(t-s)')
p('该模型称WNOJ，显式允许非零加速度，并通过q限制加速度的变化；不是强迫UAV匀速。研究背景见文献[3][4]，本项目采用R³状态特化，不声称完整实现SE(3)论文算法。')
eq(r'F(h)=\begin{bmatrix}I&hI&h^2I/2\\0&I&hI\\0&0&I\end{bmatrix}')
eq(r'Q(h)=q\begin{bmatrix}h^5/20&h^4/8&h^3/6\\h^4/8&h^3/3&h^2/2\\h^3/6&h^2/2&h\end{bmatrix}\otimes I_3')
eq(r'Q_{pp}=q\int_0^h\frac{(h-\tau)^4}{4}d\tau=\frac{qh^5}{20},\quad Q_{aa}=qh')
p('F来自对确定性ṗ=v、v̇=a积分。Q来自对白噪声积分核做协方差积分，例如位置的核是(h−τ)²/2，故其方差积分含四次幂。q单位m²/s⁵。h=0.2 s时，q=1允许单轴加速度过程变化标准差√0.2≈0.447 m/s²，q=0.01则约0.0447 m/s²。后者更强地限制机动变化，不能仅凭曲线较平就认为更可靠。')
h('7 2 两个目标克隆与全部浮标导航误差',2)
eq(r'\chi_k=[\mathbf{x}_{u,k-1}^T,\mathbf{x}_{u,k}^T,\delta\mathbf{n}_{1,k}^T,\ldots,\delta\mathbf{n}_{20,k}^T]^T')
p('稳态联合维数为2×9+20×16=338。首帧只有一个目标状态，共329维。保留前一和当前两个目标克隆，并保持全部目标与导航、不同浮标之间的交叉协方差。更老的目标克隆通过高斯边际化移出，历史观测信息仍压缩在到达因子中，不是每帧弱先验重置。')
eq(r'P_{new,rest}=FP_{last,rest},\qquad P_{new,new}=FP_{last,last}F^T+Q')
p('这两项由新状态xnew=Fxlast+w直接计算协方差得到。删除旧变量时，在协方差形式取保留变量的主子块；信息形式才对应Schur补。文献[6]支持保存跨时状态相关性的随机克隆思想，但不自动保证任意非线性边缘化实现严格一致。')
h('7 3 导航参考坐标转换',2)
eq(r'\delta\mathbf{n}_{cur}=T^{-1}(F_n\delta\mathbf{n}_{prev}-\mathbf{d}),\quad Q_{cur}=T^{-1}Q_nT^{-T}')
p('参考ESKF更新会移动坐标原点，因此把原始IMU转移与参考更新差d转换到新的误差坐标。姿态块T≈I+[dθ]×/2+[dθ]×²/12，来自旋转复合的一阶扰动展开；均值、协方差和交叉块统一变换。联合先验来自上一帧联合后验与IMU，不把参考ESKF当前P再加一次。该小角变换和高阶项截断是明确近似。')

h('八 联合MAP求解及信息保留')
h('8 1 到达先验和当前似然',2)
eq(r'J(\chi)=\frac12(\chi-\chi^-)^T(P^-)^{-1}(\chi-\chi^-)+\frac12[\mathbf{z}-h(\chi)]^TR_z^{-1}[\mathbf{z}-h(\chi)]')
p('负对数高斯密度去掉与状态无关的常数后，就是上述两个二次型。当前z包括伪距、地磁、角距以及可选FMCW距离率。上一帧原始观测已经进入后验，因此不会把整段历史测量重新加入同一先验重复计数。冻结到达先验，三次线性化迭代仅是在优化同一个后验目标。')
eq(r'\mathbf{r}_j=\mathbf{z}-h(\chi_j),\quad S_j=I+H_jP^-H_j^T,\quad K_j=P^-H_j^TS_j^{-1}')
eq(r'\chi_{j+1}=\chi^-+K_j[\mathbf{r}_j+H_j(\chi_j-\chi^-)]')
eq(r'P^+=P^--KHP^-,\qquad \Lambda^+=\Lambda^-+H^TR_z^{-1}H')
p('观测已白化，所以S中出现I。第二式可由线性化目标的正规方程和矩阵求逆引理推出；第三式是最后一次线性化的高斯协方差，只更新一次并数值对称化。算法与迭代Kalman更新或Gauss–Newton存在等价关系[4][7]，不能仅因称GP就说完全脱离KF。P只是局部高斯近似，需由NEES和覆盖率判断。')
h('8 2 原始角距与导航共同约束',2)
eq(r'H_u=[BR^T,0,0],\quad H_n=[-BR^T,0,B[\mathbf{q}_B]_\times J_r,0,0,0]')
p('B是球坐标对机体系相对位置的雅可比。来自δqB≈−Rᵀδpb+[qB]×δθ。导航变量显式进入联合状态后，雷达噪声协方差只放原始传感器噪声，不再另加Jn Pnav Jnᵀ以免重复计入。不同浮标通过共同UAV状态及联合P互相联系。')
h('8 3 薄QR和稀疏乘法为何不丢信息',2)
eq(r'G=QR,\quad \|\mathbf{r}-G\delta\|^2=\|Q^T\mathbf{r}-R\delta\|^2+\|(I-QQ^T)\mathbf{r}\|^2')
p('八条伪距的局部雅可比只有位置三列与钟差一列，故薄QR保留四个有效方向。正交补项在当前线性化下与增量无关，可以移除。这不是挑四颗好卫星。每次迭代重新线性化和压缩；每条观测只连接UAV与本浮标，所以HP和HPHᵀ只计算非零列，但P仍完整保留交叉块。最初密集版本的源码与启动日志保留，正式批次在等价压缩后全部从头运行。')
p('完整20种子无D版本与上一轮对应结果的位置和速度RMSE均值最大差异约2.1×10⁻¹¹，支持数值等价；这项复核使用完整批次，不是冒烟或门限小测试。')

h('九 未来预测与相干合成')
h('9 1 在数据真正可用后预测',2)
eq(r'h_j=0.006395+0.125+T_{compute}+\frac{j}{8}0.04,\quad j=0,\ldots,8')
p('每帧参考时刻是ADC CPI中心，采样到中心后6.395 ms才齐全。因此当前位置曲线对应有这一固定可用延迟的中心状态，不是声称在CPI中心已获得未来采样。固定处理假设取12 ms，实际时延评价取max(12 ms,实测导航均摊+ADC频谱处理+联合求解)。无D版本也处理同样ADC并承受同样等待，区别仅在是否使用速度似然。')
eq(r'\hat{\mathbf{p}}_u(h)=\hat{\mathbf{p}}_u+h\hat{\mathbf{v}}_u+\frac{h^2}{2}\hat{\mathbf{a}}_u')
eq(r'\hat{\mathbf{p}}_{b,i}(h)=\hat{\mathbf{p}}_{b,i}+h\hat{\mathbf{v}}_{b,i}+\frac{h^2}{2}[\hat R_i(\mathbf{f}_{m,i}-\hat{\mathbf{b}}_{a,i})+\mathbf{g}]')
p('UAV按当前位置、速度、加速度做二阶外推；浮标加速度从当前比力、联合修正姿态及偏置恢复。两者都不读未来观测或真值波浪相位。这轮没有运行此前的二级慢快振子预测器，也没有直接将真实慢快分量带入预测。')
h('9 2 从距离命令到复场',2)
eq(r'\hat d_i(h)=\|\hat{\mathbf{p}}_u(h)-\hat{\mathbf{p}}_{b,i}(h)\|,\quad e_{d,i}=\hat d_i-d_i')
eq(r'A_i=\sqrt{\frac{P_{tot}}{N}}\frac{\lambda_c}{4\pi d_i},\quad E=\sum_i A_i\exp(j\frac{2\pi}{\lambda_c}e_{d,i})')
eq(r'P_{rx}=|E|^2=\sum_iA_i^2+2\sum_{i<j}A_iA_j\cos[\frac{2\pi}{\lambda_c}(e_{d,i}-e_{d,j})]')
p('通信传播为单程，相位系数2π/λc；FMCW往返相位为4π/λr，两者不能混用。真实未来距离di只供评价实际通道，不反馈相位命令。20个浮标各一个各向同性等效发射点，总通信功率1 W均匀分配。载波频相同步设为理想；本轮没有功率重分配或选节点。距离与相位误差对协作波束的背景见文献[8]。')
eq(r'\eta=\frac{\sum_kP_{rx,k}}{\sum_kP_{ideal,k}},\quad G=10\log_{10}\frac{\sum_kP_{rx,k}}{\sum_k\sum_iA_{i,k}^{2}}')
p('Pideal=(ΣAi)²为同一均匀功率下理想相干功率。η是能量比，按理想功率给不同时间加权；较近帧通常权重更大。G相对非相干参考，全部20节点参与时理想等幅上限约10log10(20)=13.01 dB。高η不等价于每一帧的轨迹误差小，少数大异常可能严重拉高RMSE却只有限改变能量效率。')
h('9 3 差分几何与剩余相位误差',2)
eq(r'\delta d_i\approx\mathbf{u}_i^T(\delta\mathbf{p}_u-\delta\mathbf{p}_{b,i}),\quad \operatorname{Var}(\delta d_i-\delta d_j)=C_{ii}+C_{jj}-2C_{ij}')
eq(r'E[e^{j\Delta\phi_{ij}}]=\exp[-\frac12(\frac{2\pi}{\lambda_c})^2(C_{ii}+C_{jj}-2C_{ij})]')
p('最后一式要求零均值联合高斯距离误差，来自高斯变量特征函数，用来解释交叉协方差的重要性；它不是本轮优化发射权重的算法。当前代码保留了相对状态协方差计算，但本次五种方法的发射评价只使用点预测距离，未完成未来相位概率的严格校准或最大期望功率优化。')

h('十 完整实验设计与评价指标')
table(['方法标识','q','距离率似然'],[
('GP_noD_q1','1','不使用'),('GP_Doppler_q1','1','联合uᵀ(vu−vb)'),('GP_noD_q001','0.01','不使用'),('GP_Doppler_q001','0.01','联合uᵀ(vu−vb)'),('GP_Doppler_ignore_buoy','1','故意失配uᵀvu')],[6.8,2,7.8])
p('q是批次前固定的两档先验，不按真值运动阶段切换，也不是扫描后挑选最好值。第五种仅用于诊断忽略平台运动的后果，不是推荐模型。固定时延和含CPU时延是两种评价政策，不是两个独立随机实验。通信载频300 MHz与1 GHz也不应当成新增独立重复。')
p('种子为20261028至20261047，海况Hs=0.5与1.5 m，五类轨迹完整0至90 s，每0.2 s一帧，共451帧。共有200条物理轨迹；各轨迹运行五种方法。8个统计时段为全程、启动前5秒、启动后、15至20秒加速、35至40秒转入、40至55秒转稳、65至70秒刹停、75至90秒尾段。后七项用于解释现象，不替代全程结果。字段hover_end在直线和爬升中也只是尾段时间名。')
eq(r'20\times2\times5\times5\times2\times2\times8=32000')
p('测速指标为20×2×5×8=1600行。四个不重叠种子组各跑5种子，全部自然完成，每组约80.4分钟。原始子进程文件与最初密集版本记录均保留，没有拼接不完整前缀来冒充正式20种子研究。')
eq(r'\operatorname{RMSE}_p=\sqrt{\frac1K\sum_k\|\hat{\mathbf{p}}_k-\mathbf{p}_k\|^2}')
eq(r'J_p=\sqrt{\frac{1}{K-1}\sum_{k=2}^{K}\|\frac{\mathbf{e}_{p,k}-\mathbf{e}_{p,k-1}}{\Delta t}\|^2}')
eq(r'\operatorname{NEES}_k=\mathbf{e}_k^TP_k^{-1}\mathbf{e}_k,\quad C_{95}=\frac1K\sum_k1[\mathbf{e}_{p,k}^TP_{pp,k}^{-1}\mathbf{e}_{p,k}\le\chi^2_{3,0.95}]')
p('Jp使用相对真值的误差变化，不把真实机动直接等同抖动，但它同时包括高频波动与持续误差变化，不能单独视为纯噪声谱指标。速度误差变化率同理。六维e=[位置误差,速度误差]，正确零均值高斯协方差下NEES期望为6。覆盖率的χ²区域仅用于评价，不用于拒绝观测。')
eq(r'E_d=\sqrt{\frac{1}{KN}\sum_{k,i}(e_{d,ki}-\frac1N\sum_je_{d,kj})^2}')
p('Ed去掉同一时刻的公共距离误差，更接近相干合成需要控制的跨节点差分几何。未来UAV位置RMSE则在实际发射时刻与真值比较。全程均保留启动和异常帧，没有用后平滑美化预测轨迹。')
eq(r'\bar\Delta=\frac1{20}\sum_s(M_{new,s}-M_{old,s}),\quad CI_{half}=t_{19,0.975}\frac{s_\Delta}{\sqrt{20}}')
p('先对每个种子求RMSE或分位数，再报告20个种子的均值。95%误差分位也是各自分位数的平均，不是把所有原始样本合并后求分位。配对区间按种子差计算，不把同一状态重复的节点、帧、载频、时段行当独立样本。区间是探索性t近似，未校正多重比较，不用于声称普适最优。')

h('十一 完整结果与失效证据')
h('11 1 相同q下五轨迹的比较',2)
p('以下均为Hs=1.5 m、20个种子、完整0至90秒。无D基线是此前的联合GP，不是最早的IMM或CI模型。当前位置误差与通信载频无关。')
rows=[]
for kind,name in [('hover','悬停'),('straight','直线'),('circle','圆周'),('mixed','混合'),('climb','爬升')]:
    a=get(kind,'GP_noD_q1');b=get(kind)
    rows.append((name,f"{float(a['pos_rmse']):.4f}",f"{float(b['pos_rmse']):.4f}",f"{float(a['vel_rmse']):.4f}",f"{float(b['vel_rmse']):.4f}"))
table(['轨迹','无D位置m','有D位置m','无D速度m/s','有D速度m/s'],rows,[2.2,3.6,3.6,3.6,3.6])
p('悬停、圆周、混合三类改善明显；直线和爬升整体反而退化，不能隐藏。长直线位置RMSE配对增加0.2951±0.1484 m，爬升增加0.4869±0.1898 m，说明完整重复中也存在稳定的负面证据。')
h('11 2 混合运动的位置 速度与波动',2)
a=get('mixed','GP_noD_q1');b=get('mixed')
metrics=[('位置RMSE m','pos_rmse'),('速度RMSE m/s','vel_rmse'),('Z位置RMSE m','pZ_rmse'),('Z速度RMSE m/s','vZ_rmse'),('位置误差变化率 m/s','position_error_rate'),('速度误差变化率 m/s²','velocity_error_rate'),('未来UAV位置RMSE m','future_pos_rmse'),('去公共项距离误差RMS m','path_rmse')]
table(['指标','无D q等于1','联合D q等于1'],[(label,f"{float(a[k]):.6f}",f"{float(b[k]):.6f}") for label,k in metrics],[8.0,4.3,4.3])
p('位置、速度和误差变化率分别改善约45.5%、82.2%、65.5%，Z速度改善约75.9%。固定公共时延下，位置RMSE配对差为−0.19973±0.00718 m，速度为−0.68937±0.01464 m/s，Jp为−1.16996±0.02358 m/s。与相同q的基线比较，说明改善并非单纯收紧运动先验。')
p('下图是首种子说明，不代替20种子统计。灰色为无D q=1，红色为联合D q=1，绿色为联合D q=0.01，紫色为错误忽略浮标速度。为在Word中保留可读性，将原六面板图逐面板排放，未改变坐标或曲线数据。')
im=Image.open(DATA/'mixed_errors.png');ww,hh=im.size
for row,axis in enumerate('XYZ'):
    for col,meaning in enumerate(('位置误差','速度误差')):
        left=0 if col==0 else int(ww*.505);right=int(ww*.505) if col==0 else ww
        top=0 if row==0 else int(hh*row/3)-4;bottom=hh if row==2 else int(hh*(row+1)/3)+4
        path=WORK/f'mixed_{axis}_{col}.png';im.crop((left,top,right,bottom)).save(path)
        pic(path,f'混合运动{axis}轴{meaning} 首种子全时段')
h('11 3 置信度与两档先验的区别',2)
table(['方法或阶段','位置95%覆盖率','六维NEES'],[
('混合全程 无D q=1','96.52%','4.53'),('混合全程 联合D q=1','99.02%','2.54'),('混合全程 联合D q=0.01','98.48%','14.24'),('加速段 联合D q=1','99.80%','2.69'),('刹停段 联合D q=1','100.00%','2.36'),('加速段 联合D q=0.01','99.60%','45.19'),('刹停段 联合D q=0.01','100.00%','46.64')],[8.6,4,4])
p('位置覆盖率接近100%并不自动说明概率校准更好。q=1在混合运动中整体提示偏保守；q=0.01虽然Jp略低，但速度RMSE为0.2946 m/s，几乎是q=1的两倍，且机动六维NEES约45至47，表明速度等状态误差被低估。必须一起阅读精度、覆盖率和NEES，不能只挑最平的曲线。')
p('错误忽略浮标速度时，混合位置RMSE为1.5239 m，Z速度RMSE为0.8380 m/s，位置覆盖率仅2.37%，NEES约1053.94。此结果直接说明浮标运动项不能省略。')
h('11 4 远段测速异常与轨迹退化',2)
rows=[]
for kind,name in [('straight','直线'),('climb','爬升')]:
    aa=get(kind,'GP_noD_q1',phase='hover_end');bb=get(kind,phase='hover_end')
    rows.append((name,f"{float(aa['pos_rmse']):.3f} → {float(bb['pos_rmse']):.3f}",f"{float(aa['vel_rmse']):.3f} → {float(bb['vel_rmse']):.3f}",f"{float(bb['coverage95_position'])*100:.2f}%",f"{float(bb['nees6']):.2f}"))
table(['75至90秒','位置RMSE m','速度RMSE m/s','有D覆盖率','有D NEES'],rows,[2.6,4.2,4.2,2.8,2.8])
p('直线远段单采样ADC SNR均值约−22.73 dB，速度误差95%分位均值约0.0282 m/s，但测速RMSE均值约1.1802 m/s，呈明显长尾。少数很大的异常主导平方误差；固定σD=0.15 m/s对正常帧偏宽、对异常帧又可能过于乐观。全图选峰、远段低SNR与跳变支持检查频谱关联和观测似然，但没有独立分解噪声、几何弱方向、先验及导航耦合各自贡献。')
im=Image.open(DATA/'far_errors.png');ww,hh=im.size
for row,name in enumerate(('直线','爬升')):
    for col,meaning in enumerate(('位置误差范数','速度误差范数')):
        l=0 if col==0 else int(ww*.505);r=int(ww*.505) if col==0 else ww
        top=0 if row==0 else int(hh*.50)-4;bottom=hh if row==1 else int(hh*.50)+4
        dest=WORK/f'far_{row}_{col}.png';im.crop((l,top,r,bottom)).save(dest)
        pic(dest,f'{name}全时段{meaning} 保留晚段异常和尖峰')
h('11 5 三维轨迹和相干结果',2)
p('三维图用于展示空间路径，但各轴独立缩放，不是等比例几何；长距离路径上的小偏差应结合误差图阅读。下面提取原图的圆周、混合、爬升完整面板，未选择局部时间段。黑虚线为真值，颜色定义与前文一致。')
im=Image.open(DATA/'trajectories_3d.png');ww,hh=im.size
for index,name in [(2,'圆周'),(3,'混合'),(4,'爬升')]:
    dest=WORK/f'trajectory_{index}.png';im.crop((int(ww*index/5),int(hh*.18),int(ww*(index+1)/5),int(hh*.84))).save(dest)
    pic(dest,f'{name}三维轨迹 首种子 全部0至90秒 各轴独立缩放',width=12)
rows=[]
for freq in (300,1000):
    aa=get('mixed','GP_noD_q1',freq=freq);bb=get('mixed',freq=freq)
    rows.append((f'{freq} MHz',f"{float(aa['eta']):.6f}",f"{float(bb['eta']):.6f}",f"{float(aa['gain_db']):.3f}",f"{float(bb['gain_db']):.3f}"))
table(['通信载频','无D效率η','有D效率η','无D增益dB','有D增益dB'],rows,[3.0,3.4,3.4,3.4,3.4])
p('上述值均计入半CPI与实测部分处理成本。300 MHz混合效率提高到0.7900；1 GHz相对非相干参考提高到3.821 dB，全部20节点维持相同均匀功率，因此不是功率集中伪装成全网相干。但1 GHz效率仅0.1275，仍远离理想。剩余差分距离RMS约0.0766 m，在1 GHz对应相位量级约1.61 rad；此换算只解释敏感度，不是完整误差源分解。')
p('平均相干收益提高与远段轨迹退化可以同时发生：η是理想功率加权平均，远段权重较小；RMSE对少数大误差非常敏感；某些横向或公共误差也不完全转化为跨节点差分距离。不能用η提高反驳远段估计失效。')

h('十二 时延 适用边界与未完成工作')
table(['混合运动 Hs等于1点5米','无D q=1','联合D q=1'],[
('导航加ADC加求解均值 ms','49.57','51.81'),('发射时隙末端均值 ms','220.96','223.21'),('原帧周期 ms','200','200')],[8.6,4,4])
p('两者均超过原200 ms帧内约束。计时来自当前电脑四组并行负载，不是嵌入式最坏时间保证。ADC为20节点集中处理诊断时间，真实分布式节点并行、网络编解码、带宽和调度尚未完整建模；不能用假设的12 ms处理版本声称实时部署已经完成。')
p('物理边界包括有限谱海场与简化响应、单点机体、理想自由空间与同步、无旋翼微多普勒、海杂波、多径、扩展目标或互扰；统计边界包括固定观测方差、有限迭代、小角导航转移、高斯历史压缩与尚未校准的未来不确定度。已有论文提供机制和数学背景，不替这些本项目近似背书。')
p('本版不作为所有场景的最终默认估计器。下一步可以研究以频谱多峰和关联概率直接形成观测似然，或文献[10]的measurement t观测稳健机制，再结合角距与预测证据；同时研究弱几何方向和可执行调度。这不同于此前process t的T Trend，且本轮没有把这些未来方案冒充为第六种已实现方法。')

h('十三 源码对应与复现入口')
table(['文件或结果','实际职责'],[
('fmcw_doppler/waveform.py','复数ADC生成 Hann窗 FFT 全图峰及插值 速度测量'),
('fmcw_doppler/joint_filter.py','联合状态 相对速度因子 三次MAP QR压缩 稀疏乘法 未来预测'),
('fmcw_doppler/run_full.py','完整随机实验 五方法 两海况 五轨迹 时延与相干评分'),
('joint_gp/navigation.py','原始IMU转移 误差坐标参考及重心变换'),
('raw_sensor_maritime/physics.py','连续浮标运动 原始GNSS IMU 地磁'),
('raw_sensor_maritime/estimation.py','从原始传感器解算的导航参考'),
('imm_dpc/run_experiment.py','混合运动真值积分与插值 本轮未调用其IMM估计器'),
('gm_dpc_unified/experiment.py','悬停 直线 圆周真值 本轮未调用其旧KF作为主估计'),
('outputs/per_seed_all.csv','32000行正式合并指标'),
('outputs/summary.csv','20种子均值及区间'),
('outputs/paired_comparisons.csv','相同种子配对差及区间'),
('outputs/manifest_all.json','源hash 完整计数 物理核对与记录覆盖范围')],[7,9.6])
p('复现根目录为D:\\git\\sim2，执行fmcw_doppler/run_full.py时默认20种子；本轮为四个不重叠5种子进程，起始种子20261028、20261033、20261038、20261043，每个单BLAS线程，完整命令见同目录README。全部完成后执行report.py。重跑会更新对应输出，旧原始批次应保留备份。')
p('记录的12个执行源hash与当前文件均匹配，报告源码另有hash。两项转调轨迹文件未纳入worker起止hash清单，最终manifest补录当前hash及早于运行的修改时间，并如实披露该覆盖缺口。不能把补录说成当时已作完整加密快照。最初密集实现及日志在initial_dense和outputs/initial_dense中保留。')

h('十四 研究审查中需要回答的问题')
for question,answer in [
('为什么还称GP而不是完全脱离KF的新方法','GP提供连续时间SDE轨迹先验和信息结构；线性高斯特例与增广Kalman平滑等价。贡献在联合建模和新增真实观测机制，不在算法名称。'),
('为何不直接把FMCW速度当UAV速度','因为它是相对径向速度。浮标慢漂移和升沉都投影到视线，忽略该项的完整对照已产生明显偏差和过度自信。'),
('为何速度分辨率0点152而正常误差可以更小','FFT栅格分辨率与单目标估计精度不同，足够SNR下插值可达到亚栅格精度；零填充不增加独立信息，不能将其栅格间距当标准差。'),
('为何相干提高而直线轨迹变差','两者损失函数不同，η还按理想功率加权。必须同时报告尾段轨迹、覆盖率和实际相干，不用全程功率收益掩盖少量大误差。'),
('为什么不删掉远段的错误峰','本轮按完整样本比较，没有门控和真值限窗。数据处理的失效本身就是结果；后续需要用有依据的观测似然而非删除不利结果。'),
('这是否证明海上1 GHz协同发射已可部署','没有。1 GHz效率仍低、部分时序超200 ms、复杂海面回波与实际网络未验证，论文机制不能替本项目完成这些证明。')]:
    h(question,2);p(answer)

h('十五 正式参考文献与采用范围')
refs=[
('[1] Fossen T I, Perez T. Kalman Filtering for Positioning and Heading Control of Ships and Offshore Rigs. IEEE Control Systems Magazine, 2009, 29(6):32–46. DOI 10.1109/MCS.2009.934408.','用于低频与波频区分背景，本项目浮标响应未做船型参数移植。'),
('[2] Cohen I, Levanon N. Bistatic Radar Positioning of a Moving Target, Using Only Range and Range-Rate Measurements. IEEE Transactions on Aerospace and Electronic Systems, 2023, 59(2):1581–1589. DOI 10.1109/TAES.2022.3202758.','用于距离与距离率联合运动约束，双基地公式不能直接当单站公式。'),
('[3] Tang T Y, Yoon D J, Barfoot T D. A White-Noise-on-Jerk Motion Prior for Continuous-Time Trajectory Estimation on SE(3). IEEE Robotics and Automation Letters, 2019, 4(2):594–601. DOI 10.1109/LRA.2019.2891492.','采用R³ WNOJ特化，不借用其完整SE(3)性能结论。'),
('[4] Anderson S, Barfoot T D, Tong C H, Särkkä S. Batch nonlinear continuous-time trajectory estimation as exactly sparse Gaussian process regression. Autonomous Robots, 2015, 39:221–238. DOI 10.1007/s10514-015-9455-y.','用于SDE GP先验和线性高斯平滑等价性。'),
('[5] Kalman R E. A New Approach to Linear Filtering and Prediction Problems. Journal of Basic Engineering, 1960, 82(1):35–45. DOI 10.1115/1.3662552.','基本线性高斯递推背景，不保证本项目非线性近似自动一致。'),
('[6] Roumeliotis S I, Burdick J W. Stochastic cloning a generalized framework for processing relative state measurements. IEEE ICRA, 2002:1788–1795. DOI 10.1109/ROBOT.2002.1014801.','用于保存跨时状态相关，当前固定两克隆不等于全历史重优化。'),
('[7] Bell B M. The Iterated Kalman Smoother as a Gauss Newton Method. SIAM Journal on Optimization, 1994, 4(3):626–636. DOI 10.1137/0804035.','用于迭代更新与Gauss–Newton联系。'),
('[8] Ochiai H, Mitran P, Poor H V, Tarokh V. Collaborative Beamforming for Distributed Wireless Ad Hoc Sensor Networks. IEEE Transactions on Signal Processing, 2005, 53(11):4110–4124. DOI 10.1109/TSP.2005.857028.','用于位置及相位误差影响协同波束的背景。'),
('[9] Michalczyk J, Jung R, Weiss S. Tightly-Coupled EKF-Based Radar-Inertial Odometry. IEEE/RSJ IROS, 2022. DOI 10.1109/IROS47612.2022.9981396.','采用FMCW距离与速度同IMU紧耦合机制；原问题是静态环境中的自身里程计，非动态浮标外部跟踪。'),
('[10] Aravkin A Y, Burke J V, Pillonetto G. Robust and Trend-Following Student’s t Kalman Smoothers. SIAM Journal on Control and Optimization, 2014, 52(5):2891–2916. DOI 10.1137/130918861.','仅用于后续measurement t与process t机制区分，本轮未运行Student t版本。'),
('[11] Hoffmann F, Ritchie M, Fioranelli F, Charlish A, Griffiths H. Micro-Doppler Based Detection and Tracking of UAVs with Multistatic Radar. IEEE RadarConf, 2016. DOI 10.1109/RADAR.2016.7485236.','多站UAV观测及旋翼谱背景；是脉冲NetRAD，不是FMCW距离率融合复现。'),
('[12] Siddiq K et al. On Phase Measurement in FMCW Radar Systems. Sensor Signal Processing for Defence, 2017. DOI 10.1109/SSPD.2017.8233250.','FMCW相位、速度与无模糊性背景；本报告同时自行推导具体传播公式。')]
for citation,scope in refs:p(citation);p('采用范围 '+scope)
p('以上均指正式期刊或会议发表记录，DOI可在浏览器中以https://doi.org/加DOI打开。未把arXiv或标注Pre print的作者稿当作正式全文证据。IROS及TAES发表记录支持其方法主题，但本报告没有声称逐式复现其未完整读取的算法细节。所有本项目离散参数、保护段定义和机制组合均以源码与这里的推导为准。')

props=doc.core_properties;props.title='FMCW距离率辅助动态浮标联合跟踪UAV实验';props.subject='模型 方法 推导 完整结果与边界';props.author='项目实验说明';props.keywords='FMCW UAV 动态浮标 联合GP 相干合成'
doc.save(OUT)
print(json.dumps(dict(output=str(OUT),equations=eqcount,figures=figcount,paragraphs=len(doc.paragraphs),tables=len(doc.tables)),ensure_ascii=False))
