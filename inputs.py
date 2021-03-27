import pandas as pd
import numpy as np
from datetime import datetime
from datetime import timedelta
import os
import time as systime

## code runtime log
code_start_time = systime.perf_counter()

## Input Files

lookup_folder = r"C:\Users\####"
tariffs_folder = r"C:\Users\####"
validation_folder = r"C:\Users\####"

## tlm gsp group identifier
tlm_zone_file = lookup_folder+"\\"+"tlm_zone_prepped.xlsx"

## elexon: tlm download
tlm_data_file = lookup_folder+"\\"+"tlm.csv"

## rcrc download
rcrc_data_file = lookup_folder+"\\"+"rcrc.csv"

## bsuos download
bsuos_data_file = lookup_folder+"\\"+"bsuos_sf.csv"

## llf periods (soc data)
llf_periods_file = lookup_folder+"\\"+"llf_periods_prepped.xlsx"

## rag periods (soc data)
rag_periods_file = lookup_folder+"\\"+"rag_periods_prepped.xlsx"

## tnuos rates & times
tnuos_data_file = lookup_folder+"\\"+"tnuos_data_prepped.xlsx"

## cm rates & times
cm_data_file = lookup_folder+"\\"+"cm_prepped.xlsx"

## cfd days and rates (LCCC daily levy rates)
cfd_data_file = lookup_folder + "\\" + "cfd.xlsx"

## ia exports
ia_data_folder = lookup_folder + "\\" + "ia"
ia_data_files = []
for ia in os.listdir(ia_data_folder):
    if ia.endswith(".csv"):
        ia_data_files.append(ia)

## ia parser
ia_parser_file = lookup_folder + "\\" + "ia_parser.xlsx"

## bill folder
bill_data_folder = lookup_folder + "\\" + "hh_bills"
bill_files = []
for i in os.listdir(bill_data_folder):
        bill_files.append(bill_data_folder + "\\" + i)
        
## tariff year
year_file = "tariff_year_prepped.xlsx"

## Time Series

'''create complete time series (366 days)'''
## time series start & end day as per tariff year file
year = pd.read_excel(year_file)

## start with london time (+ 1 day before and after date range to allow full utc datetime range)
ts_start = (year.loc[0,"start_day"]-pd.Timedelta("1D")).strftime("%Y-%m-%d")
ts_end = (year.loc[0,"end_day"]+pd.Timedelta("1D")).strftime("%Y-%m-%d")
ts = pd.date_range(ts_start,ts_end,freq="30min",tz="Europe/London",name="datetime").to_frame(name="sp").sort_index()
                     
## create string for date and time (london time)
ts["date"] = ts["sp"].apply(lambda x: str(x)[:10])
ts["time"] = ts["sp"].apply(lambda x: str(x)[11:])
ts = ts.reset_index().set_index(["date","time"])

## create entries for settlement period
for day in ts.index.get_level_values("date").unique():
    counter=1
    for time in ts.loc[day,:].index.get_level_values("time").unique():
        ts.loc[(day,time),"sp"] = counter
        counter += 1

## formatting
ts.reset_index(inplace=True)
ts["sp"] = ts["sp"].astype("int64")
ts["date"] = ts["datetime"].dt.date
ts["utc_datetime"] = ts["datetime"].dt.tz_convert("utc").dt.tz_localize(None)

## "datetime" now represent utc time (tz-naive)
## BUT "date", "time", "sp" belongs to london time
ts = ts.set_index(["date","sp"], drop=True).drop("datetime",axis=1)
ts.rename(columns={"utc_datetime":"datetime"},inplace=True)

## TLM, RCRC, BSUOS

'''prep tlm'''
## get tlm zone data (based on dno)
tlm_from_dno = pd.read_excel(tlm_zone_file,sheet_name="gsp_group_from_dno",nrows=4)
tlm_from_dno["tlf_zone"] = tlm_from_dno["tlf_zone"].astype("int64")
tlm_zone = tlm_from_dno["tlf_zone"].values

## most updated tlm data from elexon*
## ie. no RF, R2, R1, SF, II data for each settlement
tlm = pd.read_csv(tlm_data_file,parse_dates=[0],dayfirst=True,index_col=[3,0,2]).drop(["Delivering","Settlement Run Type"],axis=1).sort_index(level=[0,1,2])

## filter out tlm zone for NWL's dno only
tlm = tlm.loc[(tlm_zone,),:]

## filter out 2019-20 only
tlm = tlm.loc[(slice(None),slice("2019-03-31","2020-04-02"),),:]

## swap zone name to dno name
zone_to_dno = tlm_from_dno.set_index("tlf_zone").loc[:,"dno"].to_dict()
tlm = tlm.reset_index()
tlm["Zone"] = tlm["Zone"].astype("int64").replace(zone_to_dno)
tlm.columns = ["dno","date","sp","tlm"]

## unstack tlm data of each dno into its own columns
tlm = tlm.set_index(["dno","date","sp",]).unstack(level=0)
tlm.columns = tlm.columns.map(lambda col_id: (col_id[1],col_id[0])).map("_".join)
tlm.columns.name = "tlm"

## merge dt series to settlement period
tlm = ts.merge(tlm, on=["date", "sp"], how="left")
tlm = tlm.reset_index().drop(["date","time"],axis=1).set_index("datetime")
tlm = tlm.loc["2019-04-01":"2020-03-31"]


'''form central time series'''
hh = tlm.copy()

'''prep rcrc'''
## most updated rcrc data from elexon*
## ie. no RF, R3, R2, R1, SF, II, DF for each settlement, only the latest available one
rcrc = pd.read_csv(rcrc_data_file,index_col=[0],parse_dates=True,dayfirst=True).sort_index().drop("Settlement Run Type",axis=1)

## filter out 2019-20 only
rcrc = rcrc.loc[slice("2019-03-31","2020-04-02"),:]
#### we only get 348 out of 366 days for 2019-20 rcrc

## prep and formating
rcrc.columns = ["sp","rcrc"]
rcrc.index.name = "date"
rcrc = rcrc.set_index("sp",append=True).sort_index(level=[0,1])

## merge datetime series to match to its settlement period
rcrc = ts.merge(rcrc, on=["date", "sp"], how="left")
rcrc = rcrc.reset_index().drop(["date","time"],axis=1).set_index("datetime")
rcrc = rcrc.loc["2019-04-01":"2020-03-31"]


'''form central time series'''
hh = hh.merge(rcrc["rcrc"], on = "datetime", how="outer", validate="one_to_one")

'''get bsuos'''
## note: this is SF - settlement final
bsuos = pd.read_csv(bsuos_data_file,usecols=[0,1,2],index_col=[0],parse_dates=True,dayfirst=True).sort_index()

## filter out 2019-20 only
bsuos = bsuos.loc[slice("2019-03-31","2020-04-02"),:]
#### we only get 348 out of 366 days for 2019-20 bsuos

## prep and formating
bsuos.columns = ["sp","bsuos_£/mwh"]
bsuos.index.name = "date"

bsuos["sp"] = bsuos["sp"].astype("int64")
bsuos["bsuos"] = bsuos["bsuos_£/mwh"]/10

bsuos = bsuos.set_index("sp",append=True).sort_index(level=[0,1]).drop(["bsuos_£/mwh"],axis=1)

## merge datetime series to match to its settlement period
bsuos = ts.merge(bsuos, on=["date", "sp"], how="left")
bsuos = bsuos.reset_index().drop(["date","time"],axis=1).set_index("datetime")
bsuos = bsuos.loc["2019-04-01":"2020-03-31"]


'''form central time series'''
hh = hh.merge(bsuos["bsuos"], on = "datetime", how="outer", validate="one_to_one")

## Timeseries Prep

'''create london datetime, day of week, months and time'''
## 0:MON, 1:TUE,... 6:SUN
hh["lon_time"] = hh.index.tz_localize("utc").tz_convert("europe/london")
hh["day_of_week"] = ["weekend" if (x==5 or x==6) else "weekday" for x in hh.lon_time.dt.weekday]
hh["months"] = hh.lon_time.dt.month
hh["time"] = hh.lon_time.dt.time

hh = hh.reset_index().set_index(["day_of_week","months","time"]).sort_index(level=[0,1,2])

## RAG & LLF & CM 

'''parse llf periods time data'''
for dno in ["npg","epn","lpn","spn"]:
    dno_llf = dno+"_llf"
    dno_llf = pd.read_excel(llf_periods_file,sheet_name=dno_llf)

    def time_str_to_range(time_string):
        from datetime import datetime
        from datetime import timedelta
        ## timestring format: "hh:mm – hh:mm"
        t1, t2 = time_string.split(" - ")
        ## only want time not date, set "2020-02-04" as placeholder
        t2 = (datetime.strptime("2020-02-04 "+ t2,"%Y-%m-%d %H:%M") - timedelta(minutes=1)).strftime("%H:%M")
        ##remember only comparisons need to be timezone sensitive
        return list(pd.date_range(t1, t2, freq="30min", tz="europe/london").time)

    def straighten_periods(llf_dataframe):
        import pandas as pd
        llf_dataframe["time"] = "to be filled"
        llf_dataframe["period"] = "one of four"

        llf = llf_dataframe.copy()

        for row in llf_dataframe.index:
            datacol = list(col for col in llf_dataframe.loc[[row]].dropna(axis=1).columns if col.startswith("Period"))

            for col in datacol:
                llf.loc[str(row)+str(col), :] = llf.loc[row,:]
                llf.loc[str(row)+str(col), "time"] = llf.loc[row, col]
                llf.loc[str(row)+str(col), "period"] = col

            llf.drop(row,axis=0,inplace=True)

        return llf.reset_index(drop=True).loc[:,["day_of_week","months","time","period"]]

    def flatten_time(llf_dataframe):
        import pandas as pd
        llf = llf_dataframe.copy()

        for row in llf_dataframe.index:
            counter = 1
            for time in llf_dataframe.loc[row, "time"]:
                data = llf_dataframe.loc[[row]]

                ## extract time from month range
                data["time"] = time
                data.index = [str(row)+"_"+str(counter)]

                ## append new row with sep. individual months from month range
                llf = llf.append(data)
                counter += 1

            llf.drop(row,axis=0,inplace=True)

        return llf

    def flatten_month(llf_dataframe):
        import pandas as pd
        llf = llf_dataframe.copy()

        for row in range(llf_dataframe.shape[0]):
            counter = 1
            for month in llf_dataframe.loc[row, "months"]:
                data = llf_dataframe.loc[[row]]

                ## extract months from month range
                data["months"] = month
                data.index = [str(row)+"_"+str(counter)]

                ## append new row with sep. individual month from month range
                llf = llf.append(data)
                counter += 1

            llf.drop(row,axis=0,inplace=True)

        return llf

    ## parse months
    for rows in range(dno_llf.shape[0]):
        month_range = dno_llf.loc[rows, "months"]
        if type(month_range) == str:
            dno_llf.loc[rows, "months"] = list(map(int, month_range.split("; ")))
        elif type(month_range) == int:
            dno_llf.loc[rows, "months"] = [month_range]
        else:
            print("month range data from excel is not in expected format!")

            
    ## parse time periods
    for col in list(col for col in dno_llf.columns if col.startswith("Period")):
        ## in each column, only get rows with time data
        time_rows = dno_llf.dropna(axis=0, subset=[col]).index

        for row in time_rows:
            time_data = dno_llf.loc[row, col]

            if time_data != "ALL OTHER TIMES":
                ## replace end of day from 24:00 to 23:59 (we dont match end time anyway, only start time)
                time_data = time_data.replace("24:00", "23:59")

                ## if cells contain multiple time period data
                if "\n" in time_data:
                    counter = 1
                    ## parse each time string using time_str_to_range function, then append
                    for each_time_data in time_data.split("\n"):
                        if counter == 1:
                            time_period = time_str_to_range(each_time_data)
                        else:
                            time_period.extend(time_str_to_range(each_time_data))
                        counter += 1
                    del counter
                ## if cells only contain single time period, run time_str_to_range once
                else:
                    time_period = time_str_to_range(time_data)

                ## at the end, update time string with a list of half hour start time
                dno_llf.loc[row, col] = time_period

    ## at the end of each row parsing, parse "ALL OTHER TIMES"
    for row in dno_llf.index:
        data_row = dno_llf.loc[[row],:]
        data_series = dno_llf.loc[row,:]

        if "ALL OTHER TIMES" in list(data_series.values):
            other_col = list(data_series[data_series == "ALL OTHER TIMES"].index.values)
        else:
            continue

        ## identify which column has time data
        time_col = list(col for col in data_row.dropna(axis=1).columns if col.startswith("Period"))
        for col in other_col:
            time_col.remove(col)

        ## collect all mentioned times
        mentioned_time = []
        for col in time_col:
            time_data = dno_llf.loc[row, col]
            mentioned_time.extend(time_data)

        ## "ALL OTHER TIMES" would be all remaining times outside of mentioned times
        other_time = list(pd.date_range("2020-02-04","2020-02-05",freq="30min")[:48].time)
        for mentioned in mentioned_time:
            if mentioned in other_time:
                other_time.remove(mentioned)

        ## parse "ALL OTHER TIMES" as the remaining time
        if len(other_col) == 1:
            dno_llf.at[row, other_col[0]] = other_time
        else:
            print("in this row: ", row, "more than one column with 'ALL OTHER TIMES'")

    dno_llf = straighten_periods(dno_llf)
    dno_llf = flatten_month(dno_llf)
    dno_llf = flatten_time(dno_llf)

    ## formatting 
    dno_llf = dno_llf.reset_index(drop=True).set_index(["day_of_week","months","time"]).sort_index(level=[0,1,2])
    dno_llf.columns = [dno+"_llf_period"]
    
    hh = hh.merge(dno_llf,left_index=True,right_index=True,how="left")

'''parse rag data'''
for dno in ["npg","epn","lpn","spn"]:
    dno_rag = dno+"_rag"
    dno_rag = pd.read_excel(rag_periods_file,sheet_name=dno_rag)

    def time_str_to_range(time_string):
        from datetime import datetime
        from datetime import timedelta
        ## timestring format: "hh:mm - hh:mm" *short dash not long hyphen*
        t1, t2 = time_string.split(" - ")
        t2 = (datetime.strptime("2020-02-04 "+ t2,"%Y-%m-%d %H:%M") - timedelta(minutes=1)).strftime("%H:%M")
        ##remember only comparisons need to be timezone sensitive
        return list(pd.date_range(t1, t2, freq="30min", tz="europe/london").time)

    def straighten_periods(rag_dataframe):
        import pandas as pd
        rag_dataframe["time"] = "to be filled"
        rag_dataframe["rag"] = "one of four"

        rag = rag_dataframe.copy()
        for row in rag_dataframe.index:
            datacol = list(rag_dataframe.loc[[row],:].dropna(axis=1).columns)
            datacol.remove("day_of_week")
            datacol.remove("months")
            datacol.remove("time")
            datacol.remove("rag")

            for col in datacol:
                rag.loc[str(row)+str(col), :] = rag_dataframe.loc[row,:]
                rag.loc[str(row)+str(col), "time"] = rag_dataframe.loc[row, col]
                rag.loc[str(row)+str(col), "rag"] = col
                
            rag.drop(row,axis=0,inplace=True)
            
        return rag.reset_index(drop=True).loc[:,["day_of_week","months","time","rag"]]

    def flatten_month(rag_dataframe):
        import pandas as pd
        rag = rag_dataframe.copy()
        
        for row in range(rag_dataframe.shape[0]):
            counter = 1
            ## extract months from month range
            for each_month in rag_dataframe.loc[row, "months"]:
                data = rag_dataframe.loc[[row]]
                
                ## update "month" column with "each_month" in the row's month range
                data["months"] = each_month
                data.index = [str(row)+"_"+str(counter)]
                
                ## append new row with sep. individual month from month range
                rag = rag.append(data)
                counter += 1
            del counter
            rag.drop(row,axis=0,inplace=True)
        
        return rag.reset_index(drop=True)
    
    def flatten_time(rag_dataframe):
        import pandas as pd
        rag = rag_dataframe.copy()
        for row in rag_dataframe.index:
            counter = 1
            ## extract time from month range
            for each_time in rag_dataframe.loc[row, "time"]:
                data = rag_dataframe.loc[[row]]
                
                ## updata time range with single individual "each_time"
                data["time"] = each_time
                data.index = [str(row)+"_"+str(counter)]

                ## append new row with "each_time"
                rag = rag.append(data)
                counter += 1
            del counter
            rag.drop(row,axis=0,inplace=True)
        
        return rag.reset_index(drop=True)

    ## parse months
    for rows in range(dno_rag.shape[0]):
        month_range = dno_rag.loc[rows, "months"]
        if type(month_range) == str:
            dno_rag.loc[rows, "months"] = list(map(int, month_range.split("; ")))
        elif type(month_range) == int:
            dno_rag.loc[rows, "months"] = [month_range]
        else:
            print("month range data from excel is not in expected format!")


    ## parse time periods
    for col in ["red","amber","green"]:
        ## in each column, only get rows with time data
        time_rows = dno_rag.dropna(axis=0, subset=[col]).index

        for row in time_rows:
            time_data = dno_rag.loc[row, col]

            if time_data != "ALL OTHER TIMES":
                ## replace end of day from 24:00 to 23:59 (we dont match end time anyway, only start time)
                time_data = time_data.replace("24:00", "23:59")

                ## if cells contain multiple time period data
                if "\n" in time_data:
                    counter = 1
                    ## parse each time string using time_str_to_range function, then append
                    for each_time_data in time_data.split("\n"):
                        if counter == 1:
                            time_period = time_str_to_range(each_time_data)
                        else:
                            time_period.extend(time_str_to_range(each_time_data))
                        counter += 1
                    del counter
                ## if cells only contain single time period, run time_str_to_range once
                else:
                    time_period = time_str_to_range(time_data)

                ## at the end, update time string with a list of half hour start time
                dno_rag.loc[row, col] = time_period

    ## at the end of each row parsing, parse "ALL OTHER TIMES"
    for row in dno_rag.index:
        data_row = dno_rag.loc[[row],:]
        data_series = dno_rag.loc[row,:]

        ## collect col with ALL OTHER TIMES 
        if "ALL OTHER TIMES" in list(data_series.values):
            other_col = list(data_series[data_series == "ALL OTHER TIMES"].index.values)
        else:
            continue

        ## identify which column has time data 
        ## (needs to be check! - rag periods don't have "ALL OTHER TIMES so not running this code")
        time_col = list(col for col in data_row.dropna(axis=1).columns if col in ["red","amber","green"])
        for col in other_col:
            time_col.remove(col)

        ## collect all mentioned times
        mentioned_time = []
        for col in time_col:
            time_data = dno_rag.loc[row, col]
            mentioned_time.extend(time_data)

        ## "ALL OTHER TIMES" would be all remaining times outside of mentioned times
        other_time = list(pd.date_range("2020-02-04","2020-02-05",freq="30min")[:48].time)
        for mentioned in mentioned_time:
            if mentioned in other_time:
                other_time.remove(mentioned)

        ## parse "ALL OTHER TIMES" as the remaining time
        if len(other_col) == 1:
            dno_rag.at[row, other_col[0]] = other_time
        else:
            print("in this row: ", row, ", more than one column with 'ALL OTHER TIMES'")

    dno_rag = straighten_periods(dno_rag)
    dno_rag = flatten_month(dno_rag)
    dno_rag = flatten_time(dno_rag)

    ## formatting 
    dno_rag = dno_rag.reset_index(drop=True).set_index(["day_of_week","months","time"]).sort_index(level=[0,1,2])
    dno_rag.columns = [dno+"_rag_period"]
    
    hh = hh.merge(dno_rag,left_index=True,right_index=True,how="left")

'''parse cm periods'''
cm_time = pd.read_excel(cm_data_file,sheet_name="cm_time")

def time_str_to_range(time_string):
    from datetime import datetime
    from datetime import timedelta
    ## timestring format: "hh:mm – hh:mm"
    t1, t2 = time_string.split(" - ")
    t2 = (datetime.strptime("2020-02-04 "+ t2,"%Y-%m-%d %H:%M") - timedelta(minutes=1)).strftime("%H:%M")
    ##remember only comparisons need to be timezone sensitive
    return list(pd.date_range(t1, t2, freq="30min", tz="europe/london").time)

def straighten_periods(llf_dataframe):
    import pandas as pd
    llf_dataframe["time"] = "to be filled"
    llf_dataframe["cm"] = "one of four"

    llf = llf_dataframe.copy()

    for row in llf_dataframe.index:
        datacol = list(col for col in llf_dataframe.loc[[row]].dropna(axis=1).columns if col.startswith("cm_time"))

        for col in datacol:
            llf.loc[str(row)+str(col), :] = llf.loc[row,:]
            llf.loc[str(row)+str(col), "time"] = llf.loc[row, col]
            llf.loc[str(row)+str(col), "cm"] = col

        llf.drop(row,axis=0,inplace=True)

    return llf.reset_index(drop=True).loc[:,["day_of_week","months","time","cm"]]

def flatten_time(llf_dataframe):
    import pandas as pd
    llf = llf_dataframe.copy()

    for row in llf_dataframe.index:
        counter = 1
        for time in llf_dataframe.loc[row, "time"]:
            data = llf_dataframe.loc[[row]]

            ## extract time from month range
            data["time"] = time
            data.index = [str(row)+"_"+str(counter)]

            ## append new row with sep. individual months from month range
            llf = llf.append(data)
            counter += 1

        llf.drop(row,axis=0,inplace=True)

    return llf

def flatten_month(llf_dataframe):
    import pandas as pd
    llf = llf_dataframe.copy()

    for row in range(llf_dataframe.shape[0]):
        counter = 1
        for month in llf_dataframe.loc[row, "months"]:
            data = llf_dataframe.loc[[row]]

            ## extract months from month range
            data["months"] = month
            data.index = [str(row)+"_"+str(counter)]

            ## append new row with sep. individual month from month range
            llf = llf.append(data)
            counter += 1

        llf.drop(row,axis=0,inplace=True)

    return llf

## parse months
for rows in range(cm_time.shape[0]):
    month_range = cm_time.loc[rows, "months"]
    if type(month_range) == str:
        cm_time.loc[rows, "months"] = list(map(int, month_range.split("; ")))
    elif type(month_range) == int:
        cm_time.loc[rows, "months"] = [month_range]
    else:
        print("month range data from excel is not in expected format!")


## parse time periods
for col in list(col for col in cm_time.columns if col.startswith("cm_time")):
    ## in each column, only get rows with time data
    time_rows = cm_time.dropna(axis=0, subset=[col]).index

    for row in time_rows:
        time_data = cm_time.loc[row, col]

        if time_data != "ALL OTHER TIMES":
            ## replace end of day from 24:00 to 23:59 (we dont match end time anyway, only start time)
            time_data = time_data.replace("24:00", "23:59")

            ## if cells contain multiple time period data
            if "\n" in time_data:
                counter = 1
                ## parse each time string using time_str_to_range function, then append
                for each_time_data in time_data.split("\n"):
                    if counter == 1:
                        time_period = time_str_to_range(each_time_data)
                    else:
                        time_period.extend(time_str_to_range(each_time_data))
                    counter += 1
                del counter
            ## if cells only contain single time period, run time_str_to_range once
            else:
                time_period = time_str_to_range(time_data)

            ## at the end, update time string with a list of half hour start time
            cm_time.loc[row, col] = time_period

## at the end of each row parsing, parse "ALL OTHER TIMES"
for row in cm_time.index:
    data_row = cm_time.loc[[row],:]
    data_series = cm_time.loc[row,:]

    if "ALL OTHER TIMES" in list(data_series.values):
        other_col = list(data_series[data_series == "ALL OTHER TIMES"].index.values)
    else:
        continue

    ## identify which column has time data
    time_col = list(col for col in data_row.dropna(axis=1).columns if col.startswith("cm_time"))
    for col in other_col:
        time_col.remove(col)

    ## collect all mentioned times
    mentioned_time = []
    for col in time_col:
        time_data = cm_time.loc[row, col]
        mentioned_time.extend(time_data)

    ## "ALL OTHER TIMES" would be all remaining times outside of mentioned times
    other_time = list(pd.date_range("2020-02-04","2020-02-05",freq="30min")[:48].time)
    for mentioned in mentioned_time:
        if mentioned in other_time:
            other_time.remove(mentioned)

    ## parse "ALL OTHER TIMES" as the remaining time
    if len(other_col) == 1:
        cm_time.at[row, other_col[0]] = other_time
    else:
        print("in this row: ", row, "more than one column with 'ALL OTHER TIMES'")

cm_time = straighten_periods(cm_time)
cm_time = flatten_month(cm_time)
cm_time = flatten_time(cm_time)

## formatting 
cm_time = cm_time.reset_index(drop=True).set_index(["day_of_week","months","time"]).sort_index(level=[0,1,2])
cm_time.columns = ["cm"]


'''add cm rate'''
cm_rate = pd.read_excel(cm_data_file, sheet_name="cm_provisional_rate",index_col=[0])
cm_time.replace("cm_time",cm_rate.loc["cm","p/kwh"],inplace=True)

hh = hh.merge(cm_time,left_index=True,right_index=True,how="left")

'''formatting of time series'''
hh = hh.set_index("datetime",drop=True).drop("lon_time",axis=1).sort_index().fillna(0)

## TNUOS

'''get tnuos data'''
## get tnuos costs
tnuos= pd.read_excel(tnuos_data_file,sheet_name="rates",nrows=4,usecols=[0,1])
tnuos.columns = ["dno","tnuos"]
tnuos["dno_tnuos"] = tnuos["dno"]+"_tnuos_£/kw"
tnuos.set_index("dno_tnuos",inplace=True)

## get triads time (london tz)
triads = pd.read_excel(tnuos_data_file,sheet_name="triads",nrows=3, parse_dates=True, dayfirst=True)
triads.columns = ["lon_datetime"]
triads["lon_datetime"] = triads["lon_datetime"].dt.tz_localize("europe/london")

## unpack tnuos rates to triad times
for row in triads.index:
    for dno in tnuos.index:
        triads.loc[row, dno] = tnuos.loc[dno,"tnuos"]

triads.set_index("lon_datetime",inplace=True)
triads

'''prep hh to merge'''
hh["lon_datetime"] = hh.index.tz_localize("utc").tz_convert("europe/london")
hh = hh.reset_index().set_index("lon_datetime")

hh = hh.merge(triads, on="lon_datetime", how="left").fillna(0).reset_index().set_index("datetime")
hh["lon_datetime"] = hh["lon_datetime"].dt.tz_localize(None)

## LCCC Reconciled Levy Rate

## import Reconciled Daily Levy Rate 
cfd = pd.read_excel(cfd_data_file,index_col=[0,6],parse_dates=True, dayfirst=True).sort_index(level=[0,1])
cfd = cfd.loc[(slice("2019-04-01","2020-03-31"),"Latest Settlement Run"),:].iloc[:,[0]].reset_index(level=[1],drop=True)
cfd["cfd_p/kwh"] = cfd["Reconciled Daily Levy Rate (£/MWh)"]/10
cfd = cfd.loc[:,["cfd_p/kwh"]]
cfd.index.name="date"

## form time series
days = pd.date_range("2019-04-01","2020-03-31",freq="D",name="date").to_frame(name="tempcol")
cfd = days.merge(cfd, on="date", how="left").drop("tempcol",axis=1)

## formatting
cfd.columns = ["cfd"]

## IA HH Consumption

## initialise parser
iap = pd.read_excel(ia_parser_file)
iap["mi"] = list(zip(iap.cons_type, iap.time))
iap = iap.set_index("cons")["mi"].to_dict()


'''convert and merge ia half-hourly'''
file_count = 1
for file in ia_data_files:
    file_dir = ia_data_folder+"\\"+file
    df = pd.read_csv(file_dir).drop("siteRef",axis=1).set_index(["MPAN","ConsumptionDate"])
    
    ## parse kwh & kvarh
    df.columns = pd.MultiIndex.from_tuples(df.columns.map(iap).values)
    df.columns.names = ["cons_type","time"]
    df.index.names = ["mpan","date"]
    
    ## move date from index to column then transpose (row to col, col to row)
    df = df.unstack("date").transpose().swaplevel("time","date")
    
    ## move consumption type (kwh or kvarh) to column under each mpan
    df = df.unstack("cons_type")
        
    ## parse date and time
    df.reset_index(inplace=True)
    df["datetime"] = pd.to_datetime(df['date'] + ' ' + df['time'],dayfirst=True)
    
    ## formatting
    df = df.set_index("datetime").drop(["date","time"], level=0 ,axis=1)
    df = df.transpose().reset_index()
    
    ## merging to a single main ia file
    if file_count == 1:
        ia = df.copy()
    elif file_count > 1:
        ia = ia.merge(df, on=["mpan","cons_type"], how="outer")
    
    print(file_count)
    file_count += 1
    del df

ia["mpan"] = ia["mpan"].astype("int64")


'''deal with duplicated date in the edge of ia exports'''
## get sum of duplicated edge columns
edge = ia.loc[:, (ia.columns.str[-2] == "_")].sum(axis=0).to_frame(name="sum")
edge[["dupl_datetime","dupl_id"]] = pd.DataFrame(edge.index.str.split("_").to_list(),index=edge.index)
edge.index.name="dupl_col"
edge = edge.reset_index().set_index(["dupl_datetime","dupl_id","dupl_col"]).sort_index(level=[0,1,2])

## drop the lower sum of the edge col, 
## i.e. get the most updated consumption of ia export
min_cols = []
for dt in edge.index.get_level_values(level=0).unique():
    min_cols.append(edge.loc[(dt,slice(None),slice(None)),:].idxmin().values[0][2])

ia.drop(min_cols,axis=1,inplace=True)

## remove "_x" or "_y" of the larger sum of the dupl_cols
edge = edge.reset_index(level=[0,1]).drop(min_cols)["dupl_datetime"].astype("datetime64")
ia.rename(columns = edge.to_dict(), inplace=True)
ia.set_index(["mpan","cons_type"],inplace=True)
ia.columns = pd.to_datetime(ia.columns)

del edge
del min_cols


'''set to 2019-20 datetime and set formatting'''
ia = ia.transpose().loc["2019-04-01":"2020-03-31",:]

## Bills

'''concat all 2019-20 bills'''    
for i in range(len(bill_files)):
    print(i)
    
    df = pd.read_excel(bill_files[i],
                       sheet_name="Site Summary",
                       skiprows=1,
                       usecols=["MPAN","Supply Month",
                                "CCL",
                                "MSP (kWh)","NBP (kWh)","GSP (kWh)","Consumption Charges",
                                "Management Fees","Settlement, Imbalance & Risk Charges",
                                "DUoS Rate 1","DUoS Rate 2","DUoS Rate 3",
                                "DUoS Capacity Charges","Excess Capacity Charges",
                                "DUoS Fixed Charges",
                                "Reactive Power Charges",
                                "RCRC Charges","AAHEDC Charges",
                                "BSUoS Charges","TNUoS Charges",
                                "Feed in Tariff Charges", "Contracts for Difference Charges", 
                                "Renewables Obligation Charges", 
                                "Capacity Market Charges",
                                "Daily Charges"])

    if i == 0:
        dfb = df.copy()
    elif i > 0:
        dfb = pd.concat([dfb,df])

## formatting
dfb.reset_index(drop=True, inplace=True)
dfb.columns = ["mpan","month",
               "ccl",
               "msp","gsp","nbp","commodity",
               "margin","imbalance",
               "red","amber","green",
               "capacity","excess",
               "fixed","reactive",
               "rcrc","aahedc",
               "bsuos","tnuos",
               "fit","cfd",
               "ro",
               "cm",
               "daily"
               ]

## sum all duplicated rows (reconciliation rows)
dfb = dfb.set_index(["mpan","month"]).sort_index(level=[0,1])
dfb = dfb.groupby(["mpan","month"]).sum()

## Save Prepped Files

## hh lookup
hh.to_csv(tariffs_folder + "\\" + "hh_lookup_coded.csv")

## cfd (daily lookup)
cfd.to_excel(tariffs_folder + "\\" + "cfd_coded.xlsx")

## ia (consumption)
ia.transpose().reset_index().to_csv(tariffs_folder + "\\" + "ia_coded.csv",index=False)

## bills
dfb.to_excel(validation_folder + "\\" + "bills_coded.xlsx")

## Runtime Log

code_end_time = systime.perf_counter()
print("code runtime: %0.2f minutes" %((code_end_time - code_start_time)/60))
