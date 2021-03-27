import pandas as pd
import numpy as np
import re
import time as systime
from pandas.tseries.offsets import MonthEnd
import matplotlib.pyplot as plt

## Tariffs Inputs

tariffs_folder = r"C:\Users\####"
validation_folder = r"C:\Users\####"

## dno & llfc identifier
mpan_data_file = tariffs_folder + "\\" + "mpan_data_prepped.xlsx"

## commodity data
commodity_file = tariffs_folder + "\\" + "racebank_prepped.xlsx"

## fixed rates data: imbalance, supplier margin, ro, ccl, fit, aahedc
fixed_data_file = tariffs_folder + "\\" + "fixed_rates_prepped.xlsx"

## duos rates: rag, reactive
duos_rates_file = tariffs_folder + "\\" + "duos_rates_prepped.xlsx"

## hh time series sensitive file
hh_data_file = tariffs_folder + "\\" + "hh_lookup_coded.csv"

## cfd file
cfd_data_file = tariffs_folder + "\\" + "cfd_coded.xlsx"

## llf rates
llf_rates_file = tariffs_folder + "\\" + "llf_rates_prepped.xlsx"

## ia consumption (hh of all mpans from 2019-20)
ia_file = tariffs_folder + "\\" + "ia_coded.csv"

## full mpan list
ia506_file = tariffs_folder + "\\" + "ia506.csv"

## tariff year
year_file = "tariff_year_prepped.xlsx"

'''get mpan-llfc-dno directory'''
mpan_dir = pd.read_excel(mpan_data_file,sheet_name="mpan_data",index_col=[0])
mpan_dir = mpan_dir.dropna(axis=0).astype({"llfc":"int64","dno":"str"})


'''get commodity'''
comm = pd.read_excel(commodity_file,index_col="datetime")
comm["commodity_p/kwh"] = (((comm["position_mw"]*.5-comm["racebank_mw"]*.5)*comm["n2ex"] + comm["racebank_mw"]*.5*comm["ppa"])/(comm["position_mw"]*.5))/10
comm = comm["2019-04-01":"2020-03-31"].loc[:,["commodity_p/kwh"]]
comm.columns = ["commodity"]


'''get fixed rates data'''
## imbalance, supplier margin, ro, ccl, fit, aahedc
fixed = pd.read_excel(fixed_data_file,usecols=["tariffs","p/kwh"],index_col=[0])

## duos rates: rag & reactive
duos = pd.read_excel(duos_rates_file,sheet_name="duos_rates",index_col=[0,1])
if duos.index.duplicated().any():
    print("error: 'llfc and dno pair is not unique, ie 'different dno might have the same llfc' " )

    
'''get hh data: tlm, rcrc, bsuos'''
## utc time, not localized
hhdata = pd.read_csv(hh_data_file,index_col=[0],parse_dates=True)


'''get daily data: cfd'''
## london time, not localized
cfd = pd.read_excel(cfd_data_file,index_col=[0],parse_dates=True)
cfd.index = cfd.index.date
cfd.index.name = "lon_date"


'''get llf value'''
llf = pd.read_excel(llf_rates_file).dropna(subset=["llfc"])

## parse llfc string and unpack
for row in llf.index:
    llfc_str = llf.loc[row, "llfc"]
    if type(llfc_str) == str:
        llfc_str = list(map(int, re.findall(r'\d+', llfc_str)))
        llf.at[row,"llfc"] = llfc_str
    elif type(llfc_str) == int:
        llf.at[row,"llfc"] = [llfc_str]
    
    ## unpacking
    llfc_data = llf.loc[row, "llfc"]
    llfc_count = 1
    for each_llfc in llfc_data:
        llf.loc[str(row)+"_"+str(llfc_count), :] = llf.loc[row,:]
        llf.loc[str(row)+"_"+str(llfc_count), "llfc"] = each_llfc
        llfc_count += 1
    
    llf.drop(row,axis=0,inplace=True)

llf.set_index(["dno","llfc"],inplace=True)


'''get ia'''
ia = pd.read_csv(ia_file, index_col=[0,1])
ia.columns = pd.to_datetime(ia.columns)


'''get full mpan list'''
mpan_list = pd.read_csv(ia506_file)


'''get tariff month list'''
## time series start & end day as per tariff year file
year = pd.read_excel(year_file)
ts_start = year.loc[0,"start_day"].strftime("%B %Y")
ts_end = (year.loc[0,"end_day"]-pd.Timedelta("1M")).strftime("%B %Y")

## get month list
month_range = pd.date_range(start=ts_start, end=ts_end, freq="MS").strftime("%B %Y").values

## Tariffs Formula

def tariffs(mpan, utc_start_time_string, utc_end_time_string):
    global mpan_dir, comm, fixed, duos, hhdata, cfd, llf, ia
    
    import pandas as pd
    import numpy as np
    from pandas.tseries.offsets import MonthEnd
    
    '''get dno and llfc of mpan'''
    llfc = mpan_dir.loc[mpan, "llfc"]
    dno = mpan_dir.loc[mpan, "dno"]
    
    '''initialise utc time range & create '''
    ## time range inclusive of start and end time
    t1 = utc_start_time_string
    t2 = utc_end_time_string
    
    hh = pd.date_range(t1, t2, freq="30min")
    hh = pd.DataFrame(np.random.rand(hh.shape[0]), index=hh, columns=["inicol"])
    hh.index.name = "datetime"
    hh["lon_date"] = hh.index.tz_localize("utc").tz_convert("europe/london").date
    
    '''grid multiplier'''
    ## tlm
    tlm_col = dno+"_tlm"
    hh = hh.merge(hhdata[tlm_col], on="datetime", how="left")
    
    # llf (based on llfc)
    llf_col = dno+"_llf_period"
    llf_data = hhdata.loc[:,[llf_col]].reset_index().set_index(llf_col)
    llf_rates = llf.loc[(dno, llfc),:].dropna()
    
    for period in llf_rates.index:
        llf_data.loc[period,"llf"] = llf_rates.loc[period]
    
    llf_data = llf_data.reset_index(drop=True).set_index("datetime")
    hh = hh.merge(llf_data, on="datetime", how="left")
    
    '''commodity'''
    hh = hh.merge(comm, on="datetime", how="left")
    
    '''tnuos'''
    tnuos_col = dno+"_tnuos_£/kw"
    hh = hh.merge(hhdata[tnuos_col],on="datetime",how="left")
    
    '''duos: rag'''
    ## rag rates
    rag_rates = duos.loc[(llfc, dno),["red","amber","green"]]
    
    ## rag times
    rag_col = dno+"_rag_period"
    rag = hhdata.loc[:,[rag_col]].reset_index().set_index(rag_col)
    
    ## replace rag times with rag rates (based on llfc)
    for rag_time in rag_rates.index:
        rag.loc[rag_time,"duos"] = rag_rates.loc[rag_time]
        
        if rag_time == "red":
            rag.loc[rag_time,"red"] = rag_rates.loc[rag_time]
        elif rag_time == "amber":
            rag.loc[rag_time,"amber"] = rag_rates.loc[rag_time]
    
        elif rag_time == "green":
            rag.loc[rag_time,"green"] = rag_rates.loc[rag_time]
    
    ## merge to tariff table
    rag = rag.reset_index(drop=True).set_index("datetime")
    hh = hh.merge(rag,on="datetime",how="left")
    
    '''duos: reactive'''
    hh["reactive_p/kvarh"] = duos.loc[(llfc,dno),"reactive_p/kvarh"]
    
    
    '''duos: fixed'''
    duos_fixed = duos.loc[(llfc, dno), "fixed_p/MPAN/day"]
    duos_capacity = duos.loc[(llfc, dno), "capacity_p/kVA/day"]
    duos_exceeded = duos.loc[(llfc, dno), "exceeded_p/kVA/day"]
    
    
    '''bsuos'''
    hh = hh.merge(hhdata["bsuos"], on="datetime", how="left")
    
    
    '''cm'''
    hh = hh.merge(hhdata["cm"], on="datetime", how="left")
   
    '''environmental levies'''
    ## fixed: ro, ccl, fit
    hh["ro"] = fixed.loc["ro","p/kwh"]
    hh["ccl"] = fixed.loc["ccl","p/kwh"]
    hh["fit"] = fixed.loc["fit","p/kwh"]
    
    ## daily: cfd
    hh = hh.reset_index().set_index("lon_date")
    hh = hh.merge(cfd, on="lon_date", how="left")
    hh = hh.reset_index().set_index("datetime")
    
    '''misc.'''
    ## hh: rcrc
    ## NOTE: rcrc is quoted in £/mwh
    hh = hh.merge(hhdata["rcrc"]/10, on="datetime", how="left")
    
    ## fixed: aahedc, imbalance, margin
    hh["aahedc"] = fixed.loc["aahedc","p/kwh"]
    hh["imbalance"] = fixed.loc["imbalance","p/kwh"]
    hh["margin"] = fixed.loc["margin","p/kwh"]
    
    '''total tariffs'''
    hh.drop(["inicol","lon_date"],axis=1,inplace=True)
    
    ## grid multipliers
    nbp_cols = hh[["commodity","duos","bsuos","cfd","rcrc","imbalance","margin","cm"]].sum(axis=1)
    gsp_cols = hh[["aahedc"]].sum(axis=1)
    msp_cols = hh[["ro","ccl","fit"]].sum(axis=1)
    
    nbp = hh[tlm_col]*hh["llf"]
    gsp = hh["llf"]
    
    hh["tariffs"] = nbp*nbp_cols + gsp*gsp_cols + msp_cols
    
    '''consumption'''
    hh["kwh"] = ia.loc[(mpan,"kwh"), t1 : t2]
    hh["kvarh"] = ia.loc[(mpan,"kvarh"), t1 : t2]
    hh["kw"] = hh["kwh"]/0.5
    hh["kva"] = ((hh["kwh"]**2 + hh["kvarh"]**2)**0.5)/0.5
    
    ## costs
    hh["total_cost"] = hh["tariffs"]/100 * hh["kwh"]

    return hh, tlm_col, nbp, gsp, tnuos_col, duos_fixed, duos_capacity, duos_exceeded

def month_summary(mpan, month_string):
    
    ## initialise start of month & end of month
    som = pd.Timestamp(month_string)
    eom = som + MonthEnd(1) + pd.Timedelta("1D")- pd.Timedelta("30min")
    
    ## get tariffs data
    hh, tlm_col, nbp, gsp, tnuos_col, duos_fixed, duos_capacity, duos_exceeded = tariffs(mpan, som, eom)
    
    ## prep for monthly summary
    dfc = hh.copy()
    nbp_cols = ["commodity","red","amber","green","bsuos","cfd","rcrc","imbalance","margin","cm"]
    gsp_cols = ["aahedc"]
    msp_cols = ["ro","ccl","fit"]
    
    
    '''convert unit rates to  £'''
    ## nbp, gbp and msp each have to multiply by grid multipliers
    for tariff_elements in nbp_cols:
        dfc[tariff_elements] = dfc[tariff_elements] * nbp * dfc["kwh"]/100
    for tariff_elements in gsp_cols:
        dfc[tariff_elements] = dfc[tariff_elements] * gsp * dfc["kwh"]/100
    for tariff_elements in msp_cols:
        dfc[tariff_elements] = dfc[tariff_elements] * dfc["kwh"]/100
    
    ## reactive charges
    dfc["reactive_p/kvarh"] = dfc["reactive_p/kvarh"] * dfc["kvarh"]/100
    
    ## triad charges
    dfc[tnuos_col] = dfc[tnuos_col] * dfc["kw"]
    
    ## consumption
    dfc["gsp_kwh"] = dfc["kwh"]*gsp
    dfc["nbp_kwh"] = dfc["kwh"]*nbp
    
    ## formatting
    dfc = dfc.rename(columns = {"reactive_p/kvarh":"reactive",tnuos_col:"tnuos"})
    dfc.drop([tlm_col,"llf","duos"],axis=1,inplace=True)

    
    '''monthly summary'''
    ## sum of all tariff charges
    summary = dfc.sum(axis=0)
    
    ## average hh tariffs and consumption
    summary["tariffs"] = dfc["tariffs"].mean()
    summary["kw"] = dfc["kw"].mean()
    
    ## maximum kva in a month
    summary["kva"] = dfc["kva"].max()
    
    ## formatting
    summary = summary.rename(index = {"kw":"avg_kw","kva":"max_kva","tariffs":"avg_tariffs"})
    
    # set dataframe
    summary = summary.to_frame().transpose()
    summary.index = pd.MultiIndex.from_tuples([(int(mpan),som)], names=["mpan","month"])
    
    ## add duos fixed charge (calculated as days in a month) in £
    days = summary.index.get_level_values("month").daysinmonth[0]
    summary["fixed"] = days * duos_fixed / 100

    return hh, dfc, summary

## Calculate Bills

## code start time
code_start_time = systime.perf_counter()


'''loop tariff model for each mpan'''
## initialise array for error mpans (not billed)
error_mpan = []

## start mpan count
mpan_count = 1
for mpan in mpan_list["MPAN"]:
#     print("now at mpan count:", mpan_count)
    
    month_count = 1
    for month in month_range:
        
        ## run tariff model function, return: hh rates, hh costs, and monthly summary
        try: 
            hh, dfc, summary = month_summary(mpan, month)
            
        ## if error (1. mpan not in bills, 2. some month of mpan not in bills), then: 
        except:
            print("mpan data unavailable: ",mpan)
            error_mpan.append(mpan)

        else:
            if mpan_count == 1 and month_count == 1:
                bv = summary.copy()
            else:
                bv = pd.concat([bv, summary])
        
#         print("month count: ", month_count)
        
        month_count += 1
    print("completed mpan count: ", mpan_count)
    mpan_count += 1

## formatting
bv = bv.sort_index(level=[0,1])[["ccl",
                                 "kwh","gsp_kwh","nbp_kwh","commodity",
                                 "margin","imbalance",
                                 "red","amber","green",
                                 "fixed",
                                 "max_kva","reactive",
                                 "rcrc","aahedc",
                                 "bsuos","tnuos",
                                 "fit","cfd","ro",
                                 "cm"]]
bv.rename(columns={"kwh":"msp","gsp_kwh":"gsp","nbp_kwh":"nbp"},inplace=True)

## Save Tariffs Calcs Output

## save to monthly calcs to validation folder
bv.to_excel(validation_folder + "\\" + "bv_coded.xlsx")

## Code Runtime

code_end_time = systime.perf_counter()

print("code ran for: %0.2f minutes" % ((code_end_time-code_start_time)/60))
