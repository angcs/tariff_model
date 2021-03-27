import pandas as pd
import numpy as np

## Validation Inputs

validation_folder = r"C:\Users\####"

## 2019-20 bills
bill_file = validation_folder + "\\" + "bills_coded.xlsx"

## bill calcs from tariff model
calc_file = validation_folder + "\\" + "bv_coded.xlsx"

'''get combined billed'''
bills = pd.read_excel(bill_file, index_col=[0, 1])

## create a multiindex
new_cols = []
for cols in bills.columns:
    new_cols.append(tuple(["billed",cols]))

bills.columns = pd.MultiIndex.from_tuples(new_cols).set_names(["type","charges"])

## linearise data
bills = bills.stack(level=-1)


'''get calculated bills'''
calcs = pd.read_excel(calc_file, index_col=[0, 1])

## create a multiindex
new_cols = []
for cols in calcs.columns:
    new_cols.append(tuple(["calcs",cols]))

calcs.columns = pd.MultiIndex.from_tuples(new_cols).set_names(["type","charges"])

## linearise data
calcs = calcs.stack(level=-1)

## Match Calcs & Bills

## per mpan per month
match = calcs.merge(bills, left_index=True, right_index=True, how="left")

## per charges
match_sum = match.groupby("charges").sum()
match_sum = match_sum.loc[["ccl",
                           "msp","gsp","nbp","commodity",
                           "margin","imbalance",
                           "red","amber","green",
                           "fixed",
                           "max_kva","reactive",
                           "rcrc","aahedc",
                           "bsuos","tnuos",
                           "fit","cfd","ro",
                           "cm"]]

## formatting
match = match.unstack(level=2).transpose().reorder_levels(["charges","type"]).sort_index().transpose()

## Save Matched Results

## matched files
match.to_excel(validation_folder + "\\" + "match_months.xlsx")

match_sum.to_excel(validation_folder + "\\" + "match_annual.xlsx")
