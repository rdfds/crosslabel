"""Full-centre joint connected-star workflow with explicit input stage."""
import argparse,json
from pathlib import Path
from .io import read_picasso,read_table
from .joint_workflow import STAGES,analyze_joint,format_joint_report


def main(argv=None):
    parser=argparse.ArgumentParser(prog='crosslabel-joint',description=__doc__)
    parser.add_argument('files',nargs='*',help='Picasso HDF5 files, one per label')
    parser.add_argument('--csv');parser.add_argument('--pixel-nm',type=float)
    parser.add_argument('--x',default='x');parser.add_argument('--y',default='y')
    parser.add_argument('--label-column',default='label');parser.add_argument('--scale',type=float,default=1.)
    parser.add_argument('--input-stage',choices=STAGES,required=True)
    parser.add_argument('--bootstrap',type=int,default=999);parser.add_argument('--seed',type=int,default=0)
    parser.add_argument('--radius',type=float,default=30.);parser.add_argument('--json')
    args=parser.parse_args(argv)
    if args.json and Path(args.json).exists():raise FileExistsError('output exists; refusing to overwrite evidence')
    if bool(args.csv)==bool(args.files):parser.error('provide either --csv or per-label HDF5 files')
    if args.csv:xy,labels,names=read_table(args.csv,args.x,args.y,args.label_column,args.scale)
    else:
        if args.pixel_nm is None:parser.error('--pixel-nm is required for HDF5 input')
        xy,labels=read_picasso(args.files,args.pixel_nm);names=args.files
    result=analyze_joint(xy,labels,input_stage=args.input_stage,seed=args.seed,B=args.bootstrap,radius=args.radius)
    result['label_names']=names
    print(format_joint_report(result))
    if args.json:
        with open(args.json,'x') as handle:json.dump(result,handle,indent=2,allow_nan=False)
    return 0 if result['status']=='estimated' else 2


if __name__=='__main__':raise SystemExit(main())
