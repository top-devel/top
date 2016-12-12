source $HOME/TOP_2016/top/BUILD/activate-top.sh
#top-build --parser cesam1D.eq --model=cesam --order=order-cesam-1D.inc         
#top-build --parser cesam1D_lag_l0.eq --model=cesam --order=order-cesam-1D_l0.inc
top-build --parser cesam1D_eul_l0.eq --model=cesam --order=order-cesam-1D_l0.inc
top-build --parser cesam1D_eul_iso_l0.eq --model=cesam --order=order-cesam-1D_l0.inc
#top-build --parser cesam1D_lag_iso_l0.eq --model=cesam --order=order-cesam-1D_l0.inc
python ./run_model.py
