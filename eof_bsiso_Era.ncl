;load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/home/users/bi913747/NCL/China_map/script/cnmask.ncl"
load "/home/users/bi913747/NCL/SouthChinaSea/scsstation.ncl"
begin
;读取850hpa纬向风、OLR、降水资料，时间是1982-2011
;it_s=19810101
it_s=20010101
it_e=20101231

   yrStrt  = it_s/10000
   yrLast  = it_e/10000
;;;;;;;;;;;;;data eof;;;;;;;;;;;;;;
    f_olr =addfile("/gws/nopw/j04/klingaman/bi913747/OBS/ERAinterim/TTR/ERA-Interim.daily.1979-2016.TTR.nc","r")
    time_olr=f_olr->t
    YYYYMMDD_olr=cd_calendar(time_olr,-2)      ;这里用作多变量分析的数据可以在读取时限制到需要分析的范围，
    numb = dimsizes(YYYYMMDD_olr)
    rec_s=ind(it_s.eq.YYYYMMDD_olr)            ;因为只是对时间维操作，分析范围外的数据读取只会增加计算量
    rec_e=ind(it_e.eq.YYYYMMDD_olr)
;    olr=short2flt(f_olr->TTR(rec_s:rec_e,0,{-10:40},{40:160}))
    olr=f_olr->TTR(rec_s:rec_e,0,{-10:40},{40:160})
    delete(rec_s)
    delete(rec_e)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    f_uwnd =addfile("/gws/nopw/j04/klingaman/bi913747/OBS/ERAinterim/U/ERA-Interim.daily.1979-2016.U850.nc","r")
    time_uwnd=f_uwnd->t
    YYYYMMDD_uwnd=cd_calendar(time_uwnd,-2)      ;这里用作多变量分析的数据可以在读取时限制到需要分析的范围，
    numb = dimsizes(YYYYMMDD_uwnd)
    rec_s=ind(it_s.eq.YYYYMMDD_uwnd)            ;因为只是对时间维操作，分析范围外的数据读取只会增加计算量
    rec_e=ind(it_e.eq.YYYYMMDD_uwnd)
;    uwnd=short2flt(f_uwnd->U(rec_s:rec_e,0,{-10:40},{40:160}))
    uwnd=f_uwnd->U(rec_s:rec_e,0,{-10:40},{40:160})
    delete(rec_s)
    delete(rec_e)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    f_vwnd =addfile("/gws/nopw/j04/klingaman/bi913747/OBS/ERAinterim/V/ERA-Interim.daily.1979-2016.V850.nc","r")
    time_vwnd=f_vwnd->t
    YYYYMMDD_vwnd=cd_calendar(time_vwnd,-2)      ;这里用作多变量分析的数据可以在读取时限制到需要分析的范围，
    numb = dimsizes(YYYYMMDD_vwnd)
    rec_s=ind(it_s.eq.YYYYMMDD_vwnd)            ;因为只是对时间维操作，分析范围外的数据读取只会增加计算量
    rec_e=ind(it_e.eq.YYYYMMDD_vwnd)
;    vwnd=short2flt(f_vwnd->V(rec_s:rec_e,0,{-10:40},{40:160}))
    vwnd=f_vwnd->V(rec_s:rec_e,0,{-10:40},{40:160})
    time=f_vwnd->t(rec_s:rec_e)
    delete(rec_s)
    delete(rec_e)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;得到anomaly data ，先remove annual cycle,然后remove interannual variability
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   nhar=3
      TIME    = cd_calendar(time, 0)             ; type float 
      year    = floattointeger( TIME(:,0) )          ;依次取出时间变量的年、月、日
      month   = floattointeger( TIME(:,1) ) 
      day     = floattointeger( TIME(:,2) ) 
      ddd     = day_of_year(year, month, day)    ; ddd中每个值代表对应日期在对应年份中的日期顺序
      yyyyddd = year*1000 + ddd 

      olr_ClmDay = clmDayTLL(olr, yyyyddd)     ; daily climatology at each grid point
      olr_ClmDay_sm = smthClmDayTLL(olr_ClmDay, nhar)
      uwnd_ClmDay = clmDayTLL(uwnd, yyyyddd)     ; daily climatology at each grid point
      uwnd_ClmDay_sm = smthClmDayTLL(uwnd_ClmDay, nhar)
      vwnd_ClmDay = clmDayTLL(vwnd, yyyyddd)     ; daily climatology at each grid point
      vwnd_ClmDay_sm = smthClmDayTLL(vwnd_ClmDay, nhar)
   
      olr_Anom_sm   = calcDayAnomTLL (olr, yyyyddd, olr_ClmDay_sm)
      uwnd_Anom_sm  = calcDayAnomTLL (uwnd, yyyyddd, uwnd_ClmDay_sm)
      vwnd_Anom_sm  = calcDayAnomTLL (vwnd, yyyyddd, vwnd_ClmDay_sm)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   
;   ;下面remove interannual variability(running mean of the preceeding 120 days)
;   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  dims_olr_anom=dimsizes(olr_Anom_sm)
   
   mt=dims_olr_anom(0)
   mlat=dims_olr_anom(1)
   mlon=dims_olr_anom(2)
   ave=new((/mlat,mlon/),float)
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;;;
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   olr_ano=new((/mt,mlat,mlon/),float)
   do i=0,118
      olr_ano(i,:,:)=olr_Anom_sm(i,:,:)
   end do
   do i=119,mt-1
      ave=0
      do j=0,119
        ave=ave+olr_Anom_sm(i-j,:,:)
     end do
     ave=ave/120.0
   olr_ano(i,:,:)=olr_Anom_sm(i,:,:)-ave
   end do
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;;;
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   uwnd_ano=new((/mt,mlat,mlon/),float)
   do i=0,118
      uwnd_ano(i,:,:)=uwnd_Anom_sm(i,:,:)
   end do
   do i=119,mt-1
      ave=0
      do j=0,119
        ave=ave+uwnd_Anom_sm(i-j,:,:)
     end do
     ave=ave/120.0
   uwnd_ano(i,:,:)=uwnd_Anom_sm(i,:,:)-ave
   end do
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   vwnd_ano=new((/mt,mlat,mlon/),float)
   do i=0,118
      vwnd_ano(i,:,:)=vwnd_Anom_sm(i,:,:)
   end do
   do i=119,mt-1
      ave=0
      do j=0,119
        ave=ave+vwnd_Anom_sm(i-j,:,:)
     end do
     ave=ave/120.0
   vwnd_ano(i,:,:)=vwnd_Anom_sm(i,:,:)-ave
   end do
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
   copy_VarMeta(olr, olr_ano)
   copy_VarMeta(uwnd,uwnd_ano)   
   copy_VarMeta(vwnd,vwnd_ano)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   
   ;提取夏季数据
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   date = cd_calendar(time&t,0)                          ;date变量用来存储19820101到20111231的时间 
   month  = tointeger(date(:,1))
   nt=ind(month.ge.5.and.month.le.10)                       ;选择5月到10月的数据
   olr_sum=olr_ano(nt,:,:)
   uwnd_sum=uwnd_ano(nt,:,:)
   vwnd_sum=vwnd_ano(nt,:,:)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;对夏季 太平洋区域数据 进行标准化处理，得到olr_WP_norma,uwnd_WP_norma
   ;先限制OLR，U850在太平洋区域，然后每个格点标准差，再求区域平均标准差，变量场除以平均标准差即可。
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      latS_WP    = -10
      latN_WP    =  40
      lonW_WP    =  40
      lonE_WP    =  160
   
      olr_WP=olr_sum(:,{latS_WP:latN_WP},{lonW_WP:lonE_WP})
      uwnd_WP=uwnd_sum(:,{latS_WP:latN_WP},{lonW_WP:lonE_WP})        ;提取太平洋的OLR和纬向风数据
      vwnd_WP=vwnd_sum(:,{latS_WP:latN_WP},{lonW_WP:lonE_WP})
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
      dims=dimsizes(olr_WP)      
      std_field=new((/3,dims(1),dims(2)/),float)                                   ;求标准差场
      std_field(0,:,:)=dim_stddev_n_Wrap(olr_WP,0)
      std_field(1,:,:)=dim_stddev_n_Wrap(uwnd_WP,0)
      std_field(2,:,:)=dim_stddev_n_Wrap(vwnd_WP,0)
   
      a=dim_avg_n_Wrap(dim_avg_n_Wrap(std_field(0,:,:), 0), 0)
      b=dim_avg_n_Wrap(dim_avg_n_Wrap(std_field(1,:,:), 0), 0)
      c=dim_avg_n_Wrap(dim_avg_n_Wrap(std_field(2,:,:), 0), 0)
   
print(a)
print(b)
print(c)

       olr_WP_norma=olr_WP/a
       uwnd_WP_norma=uwnd_WP/b
       vwnd_WP_norma=vwnd_WP/c


exit
end 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   
;   ;对太平洋MV-EOF分析，提取PC1和PC2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      neof=4
;   
;      dim_WP_norma   = dimsizes(olr_WP_norma)                                     
;   
;   cdata   = new ( (/2*dim_WP_norma(1),2*dim_WP_norma(2),dim_WP_norma(0)/), float)   ;新建一个数组，存放OLR和uwnd数据，三个维度分别是：lat、lon、time
;      do m1=0,dim_WP_norma(1)-1
;       do m2=0,dim_WP_norma(2)-1
;           cdata(m1,       m2,       :)=(/olr_WP_norma(:,m1,m2)/)                  
;           cdata(m1+dim_WP_norma(1),  m2+dim_WP_norma(2),  :)=(/uwnd_WP_norma(:,m1,m2)/)      
;       end do 
;      end do
;   
;   eof_cdata    = eofunc(cdata   , neof, False)    ;对cdata经验正交分解，返回三维数组，对应的是不同特征值相应特征向量的分布
;   eof_ts_cdata = eofunc_ts(cdata,eof_cdata,False) ;计算前两个模态对应的时间序列
;   printVarSummary(eof_cdata)
;   printVarSummary(eof_ts_cdata)
;   
;   PC1           =   eof_ts_cdata(2,:)
;   PC2           =  -eof_ts_cdata(3,:)
;   
;   olr_eof   =new((/2,dim_WP_norma(1),dim_WP_norma(2)/),float)
;   ;olr_eof2   =new((/dim_WP_norma(1),dim_WP_norma(2)/),float)
;   uwnd_eof  =new((/2,dim_WP_norma(1),dim_WP_norma(2)/),float)
;   ;uwnd_eof2  =new((/dim_WP_norma(1),dim_WP_norma(2)/),float)
;   vwnd_eof  =new((/2,dim_WP_norma(1),dim_WP_norma(2)/),float)
;   ;vwnd_eof2  =new((/dim_WP_norma(1),dim_WP_norma(2)/),float)
;   
;   olr_eof(0,:,:)  =  eof_cdata(2,0:dim_WP_norma(1)-1,0:dim_WP_norma(2)-1)   ;提取第一、第二模态的OLR场
;   olr_eof(1,:,:)  = -eof_cdata(3,0:dim_WP_norma(1)-1,0:dim_WP_norma(2)-1)
;   
;   uwnd_eof(0,:,:) =  eof_cdata(2,dim_WP_norma(1):2*dim_WP_norma(1)-1,dim_WP_norma(2):2*dim_WP_norma(2)-1) ;提取第一、第二模态纬向风场
;   uwnd_eof(1,:,:) = -eof_cdata(3,dim_WP_norma(1):2*dim_WP_norma(1)-1,dim_WP_norma(2):2*dim_WP_norma(2)-1)
;   
;   
;   vwnd_eof(0,:,:) = regCoef_n(PC1,vwnd_WP_norma, 0, 0)          ;vwnd_WP_norma与PC1、PC2序列做回归，得到vwnd_eof1,vwnd_eof2
;   vwnd_eof(1,:,:) = regCoef_n(PC2,vwnd_WP_norma, 0, 0)
;   
;   ;对上面各变量经度、纬度赋予属性
;   dim_olr_eof   = dimsizes(olr_eof)
;   mlat          = dim_olr_eof(1)
;   mlon          = dim_olr_eof(2)
;   ceof          = new((/3,2,mlat,mlon/),float)
;   do n=0,1
;        ceof(0,n,:,:) = olr_eof(n,:,:)              ; precip
;        ceof(1,n,:,:) = uwnd_eof(n,:,:)             ; u
;        ceof(2,n,:,:) = vwnd_eof(n,:,:)             ; v
;   end do
;   
;   ceof!0   = "var"
;   ceof!1   = "eof"
;   ceof!2   = "lat"
;   ceof!3   = "lon" 
;   
;   lat      =olr_WP&lat
;   lon      =olr_WP&lon
;   lat@units="degrees north"
;   lon@units="degrees east"
;   
;   ceof&lat = lat
;   ceof&lon = lon
;   printVarSummary(ceof) 
;   printMinMax(ceof, 0)
;   
;   ;plot
;   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;   latS = -10
;   latN = 40
;   lonW = 40
;   lonE = 160
;   
;   
;    wks  = gsn_open_wks("pdf","eof_bsiso")
;    plot = new(2,graphic)              ; create graphic array
;    res                      = True 
;    
;    ;一页多图与多图层设置       
;     res@gsnDraw              = False        ; don't draw yet
;     res@gsnFrame             = False        ; don't advance frame yet
;   
;   
;   ;地图参数设置
;     res@mpFillOn             = False        ; 地图不填充
;     res@mpMinLatF            = latS         ; 地图最小纬度
;     res@mpMaxLatF            = latN         ; 地图最大纬度
;     res@mpMinLonF            = lonW
;     res@mpMaxLonF            = lonE
;     res@gsnAddCyclic         =False
;     ;等值线设置
;     res@cnFillOn             = True                    ; 打开填色图
;     res@cnFillPalette        = "BlueWhiteOrangeRed"    ; 设置填充颜色
;     res@cnLinesOn            = False                   ; 不绘制等值线
;     res@cnLineLabelsOn       = False                   ; 不打开等值线标签
;   
;     res@lbLabelBarOn         = False              ; 关闭每个图形各自的色标，以公用一个色标
;   
;     res@gsnScalarContour     = True                ; 
;     res@gsnMajorLatSpacing   = 15                  ;主纬度间隔是15
;     res@gsnMajorLonSpacing   = 30                  ;主经度间隔是60
;   
;   
;   
;   
;   ;mnmxint = nice_mnmxintvl( min(x) , max(x), 16, False),设定等值线之前，可以先判断最小最大值，先用该函数判断，再重新设定
;   ;设定等值线最大最小值，及间隔
;     res@cnLevelSelectionMode = "ExplicitLevels"                     ;设置绘制指定等值线
;     ; res@cnMinLevelValF       =  -0.09         ; -100; mnmxint(0)    ;
;     ; res@cnMaxLevelValF       =   0.07       ;   80; mnmxint(1)
;     ; res@cnLevelSpacingF      =   0.01         ;   20; mnmxint(2)
;   
;     res@cnLevels             =(/-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0.01,0.02,0.03,0.04,0.05,0.06,0.07/)
;     res@cnFillColors         =(/10,20,25,30,40,50,60,70,80,125,185,195,205,215,225,235,245/)
;   
;   
;   res@tmXBOn                     =True
;   res@tmXBLabelDeltaF            =-0.6
;   res@tmXBLabelFontHeightF       =0.005
;   res@tmXBLabelFontThicknessF    =0.005
;   
;   res@tmXBMajorLengthF           =0.003
;   res@tmXBMajorOutwardLengthF    =0.003
;   
;   res@tmXBMinorOn                =True
;   res@tmXBMinorLengthF           =0.0
;   res@tmXBMinorOutwardLengthF    =0.0
;   
;   res@tmYLOn                     =True
;   res@tmYLLabelDeltaF            =-0.6
;   res@tmYLLabelFontHeightF       =0.005
;   res@tmYLLabelFontThicknessF    =0.005
;   
;   res@tmYLMajorLengthF           =0.003
;   res@tmYLMajorOutwardLengthF    =0.003
;   
;   res@tmYLMinorOn                =True
;   res@tmYLMinorLengthF           =0.0
;   res@tmYLMinorOutwardLengthF    =0.0
;   
;   
;   res@tmXTOn                    =False
;   res@tmYROn                    =False
;   
;   
;   res@vpWidthF                   =0.15
;   res@vpHeightF                  =0.15
;   ;res@gsnPanelXWhiteSpacePercent = 1.2
;   
;   res@lbBoxEndCapStyle           = "TriangleBothEnds"
;   
;   ;矢量箭头设置
;     res@vcMinDistanceF            = 0.01            ; 控制vector的稀疏
;     res@vcRefMagnitudeF           = 0.05            ; 定义矢量参考箭头的量
;     res@vcRefLengthF              = 0.025           ; 定义矢量参考箭头的长度
;     res@vcRefAnnoOrthogonalPosF   = -1.0            ; 参考箭头的垂直位置
;     res@vcRefAnnoArrowLineColor   = "black"         ; 设置矢量箭头的颜色
;     res@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref
;   
;   
;    resP                     = True         ; modify the panel plot
;    resP@gsnPanelLabelBar    = True 
;    resP@lbLabelFontHeightF  = 0.0095
;    resP@pmLabelBarWidthF    = 0.4          ; 设置色标的宽度
;    resP@pmLabelBarHeightF   = 0.05         ; 设置色标的高度
;    ;presP@gsnPanelFigureStringsFontHeightF = 0.0125
;    ; resP@gsnPanelDebug       =True
;    ; resP@gsnPanelXF          =(/0.0629,0.0629/)   
;   
;   
;   
;    ; resP@gsnPanelRowSpec             = True
;   
;   
;   do i=0,1
;     if (i.eq.0) then
;       res@gsnLeftString             ="a) EOF1 18.4%"
;       res@gsnLeftStringFontHeightF  =0.005
;       plot(i) = gsn_csm_vector_scalar_map(wks,ceof(1,0,:,:),ceof(2,0,:,:),ceof(0,0,:,:),res)
;      else 
;       res@gsnLeftString             ="b) EOF2 11.6%"
;       res@gsnLeftStringFontHeightF  =0.005
;       plot(i) = gsn_csm_vector_scalar_map(wks,ceof(1,1,:,:),ceof(2,1,:,:),ceof(0,1,:,:),res)
;      end if 
;   end do
;   
;   ; plot(2) = gsn_csm_xy(wks,PC&time,PC(mode|:,time|:),resxy)
;   
;   ;resP@gsnPanelFigureStrings= "WPISO"
;   resP@gsnPanelMainString  = yrStrt+"-"+yrLast
;   resP@gsnPanelMainFontHeightF = 0.02
;   gsn_panel(wks,plot,(/1,2/),resP)
;   frame(wks)
;   
;   end
;   
;   
