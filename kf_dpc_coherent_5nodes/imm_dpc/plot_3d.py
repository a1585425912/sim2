"""3D estimated trajectories and error trajectories from saved first seed."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

OUT=Path(__file__).resolve().parent/'outputs'
KINDS={'hover':'悬停','straight':'直线','circle':'圆周','mixed':'混合运动'}
METHODS={'CV_original':('CV+DPC','#3b82b8'), 'CT':('CT+DPC','#db8c28'),
         'IMM_SC':('双模型IMM+DPC','#339477'), 'IMM':('三模型IMM+DPC','#c84c65')}
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],
                     'axes.unicode_minus':False,'font.size':10})


def main():
    d=np.load(OUT/'first_seed.npz');mask=d['time_s']>=5
    for errors in (False,True):
        fig=plt.figure(figsize=(15,12))
        for index,(kind,label) in enumerate(KINDS.items()):
            ax=fig.add_subplot(2,2,index+1,projection='3d')
            truth=d[f'{kind}__truth_p'][mask]
            for method,(name,color) in METHODS.items():
                estimate=d[f'{kind}__{method}__state'][mask,0,:3]
                p=estimate-truth if errors else estimate
                ax.plot(*p.T,color=color,lw=1.2,alpha=.85)
            if errors:
                ax.scatter([0],[0],[0],marker='*',s=120,color='black',depthshade=False)
            else:
                ax.plot(*truth.T,color='black',ls='--',lw=2.2)
                ax.scatter(*truth[0],color='black',marker='o',s=35,depthshade=False)
                ax.scatter(*truth[-1],color='black',marker='*',s=90,depthshade=False)
            prefix='误差 ' if errors else ''
            ax.set_xlabel(prefix+'X (m)',labelpad=8)
            ax.set_ylabel(prefix+'Y (m)',labelpad=8)
            ax.set_zlabel(prefix+'Z (m)',labelpad=9)
            ax.set_title(label,pad=12)
            ax.view_init(elev=24,azim=-58)
            ax.set_box_aspect((1.1,1,.65))
            for axis in (ax.xaxis,ax.yaxis,ax.zaxis):
                axis.set_major_locator(MaxNLocator(5));axis.get_major_formatter().set_useOffset(False)
        handles=[Line2D([0],[0],color='black',ls='--',lw=2,label='零误差位置 ★' if errors else '真值（●起点 / ★终点）')]
        handles += [Line2D([0],[0],color=c,lw=2,label=n) for n,c in METHODS.values()]
        fig.legend(handles=handles,loc='lower center',ncol=3,frameon=False,bbox_to_anchor=(.5,.015))
        title='3D误差轨迹：估计坐标减真值' if errors else '3D运动轨迹：各模型估计与真值对照'
        fig.suptitle(title+'\n首种子20261028，节点1，5–90s；各坐标轴独立缩放，非等比例几何',fontsize=15,y=.98)
        fig.subplots_adjust(left=.03,right=.93,top=.90,bottom=.1,wspace=.12,hspace=.20)
        name='trajectory_errors_3d' if errors else 'trajectory_fits_3d'
        fig.savefig(OUT/f'{name}.png',dpi=170)
        fig.savefig(OUT/f'{name}.svg')
        plt.close(fig)
    print('Saved 3D trajectory and error plots, PNG + SVG.')

if __name__=='__main__':main()
