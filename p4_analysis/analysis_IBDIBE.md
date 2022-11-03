IBD/IBE Analysis
================

``` r
library("here")
library("sjPlot")
library("tidyverse")
library("lme4")
library("viridis")
library("lmerTest")
library("ggplot2")
library("gridExtra")
library("gt")
library("ggthemes")

source(here("p4_analysis", "analysis_functions.R"))

p4path <- here("p4_analysis", "outputs")
p3path <- here("p3_methods", "outputs")
```

# 1. MMRR

## 1.1 Individual sampling

### 1.1.1 Summary plots

``` r
mmrr_ind <- format_mmrr(here(p3path, "mmrr_indsampling_results.csv"))

# overall error
summary_hplot(mmrr_ind, "ratio_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
summary_hplot(mmrr_ind, "env_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
summary_hplot(mmrr_ind, "geo_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

``` r
# bias
summary_hplot(mmrr_ind, "ratio_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-2-4.png)<!-- -->

``` r
summary_hplot(mmrr_ind, "comboenv_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-2-5.png)<!-- -->

``` r
summary_hplot(mmrr_ind, "geo_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-2-6.png)<!-- -->

### 1.1.2 Model summaries

``` r
run_lmer(mmrr_ind, "ratio_ae", filepath = here(p4path, "MMRR_individual_ratioerr.csv"))
```

    ## <div id="fixjlqekcp" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #fixjlqekcp .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #fixjlqekcp .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #fixjlqekcp .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #fixjlqekcp .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #fixjlqekcp .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #fixjlqekcp .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #fixjlqekcp .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #fixjlqekcp .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #fixjlqekcp .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Predictors</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Fixed Effects</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Sum Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Mean Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">NumDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">DenDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">F value</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Pr(&gt;F)</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">nsamp</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;0.0003</td>
    ## <td class="gt_row gt_right">8.0000</td>
    ## <td class="gt_row gt_right">8.0000</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">1531.810409</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.9 &times; 10<sup class='gt_super'>&minus;319</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.1078</td>
    ## <td class="gt_row gt_right">1.5728</td>
    ## <td class="gt_row gt_right">0.5243</td>
    ## <td class="gt_row gt_right">3</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">100.385844</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.4 &times; 10<sup class='gt_super'>&minus;64</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #F0F5EF; color: #000000;">0.0126</td>
    ## <td class="gt_row gt_right">0.6090</td>
    ## <td class="gt_row gt_right">0.6090</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">116.601181</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">4.4 &times; 10<sup class='gt_super'>&minus;27</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #DAF0D4; color: #000000;">0.0523</td>
    ## <td class="gt_row gt_right">10.5133</td>
    ## <td class="gt_row gt_right">10.5133</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">2013.030919</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F6F6; color: #000000;">&minus;0.0020</td>
    ## <td class="gt_row gt_right">0.0155</td>
    ## <td class="gt_row gt_right">0.0155</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">2.971532</td>
    ## <td class="gt_row gt_right">0.085</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">H</td>
    ## <td class="gt_row gt_right" style="background-color: #F2F6F1; color: #000000;">0.0096</td>
    ## <td class="gt_row gt_right">0.3544</td>
    ## <td class="gt_row gt_right">0.3544</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">67.857541</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.9 &times; 10<sup class='gt_super'>&minus;16</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F4F1F4; color: #000000;">&minus;0.0090</td>
    ## <td class="gt_row gt_right">0.3123</td>
    ## <td class="gt_row gt_right">0.3123</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">59.800103</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.1 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="8"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'pbkrtest.limit = 15360' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 15360)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 15360' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 15360)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## <div id="tpzhtafrqn" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #tpzhtafrqn .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #tpzhtafrqn .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #tpzhtafrqn .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #tpzhtafrqn .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #tpzhtafrqn .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Contrast</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Estimate</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">SE</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Z ratio</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">p</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - grid<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #F1EAF1; color: #000000;">&minus;0.0046</td>
    ## <td class="gt_row gt_right">0.0016</td>
    ## <td class="gt_row gt_right">&minus;2.7960</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.02658997<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">envgeo - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F5F6; color: #000000;">&minus;0.0009</td>
    ## <td class="gt_row gt_right">0.0016</td>
    ## <td class="gt_row gt_right">&minus;0.5319</td>
    ## <td class="gt_row gt_right">0.95132229</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - trans<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #AF8DC3; color: #000000;">&minus;0.0249</td>
    ## <td class="gt_row gt_right">0.0016</td>
    ## <td class="gt_row gt_right">&minus;15.0697</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">grid - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #EEF5EC; color: #000000;">0.0037</td>
    ## <td class="gt_row gt_right">0.0016</td>
    ## <td class="gt_row gt_right">2.2641</td>
    ## <td class="gt_row gt_right">0.10653142</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">grid - trans<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #C4A7D1; color: #000000;">&minus;0.0202</td>
    ## <td class="gt_row gt_right">0.0016</td>
    ## <td class="gt_row gt_right">&minus;12.2737</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">rand - trans<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #B392C6; color: #000000;">&minus;0.0240</td>
    ## <td class="gt_row gt_right">0.0016</td>
    ## <td class="gt_row gt_right">&minus;14.5378</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.05</td>
    ##     </tr>
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">**</sup> p &lt; 0.001</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

``` r
run_lmer(mmrr_ind, "geo_ae", filepath = here(p4path, "MMRR_individual_geoerr.csv"))
```

    ## <div id="fixjlqekcp" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #fixjlqekcp .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #fixjlqekcp .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #fixjlqekcp .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #fixjlqekcp .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #fixjlqekcp .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #fixjlqekcp .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #fixjlqekcp .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #fixjlqekcp .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #fixjlqekcp .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Predictors</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Fixed Effects</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Sum Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Mean Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">NumDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">DenDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">F value</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Pr(&gt;F)</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">nsamp</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;0.0001</td>
    ## <td class="gt_row gt_right">1.2007</td>
    ## <td class="gt_row gt_right">1.2007</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">1.422367e+03</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">9.2 &times; 10<sup class='gt_super'>&minus;298</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.0525</td>
    ## <td class="gt_row gt_right">1.3073</td>
    ## <td class="gt_row gt_right">0.4358</td>
    ## <td class="gt_row gt_right">3</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">5.161977e+02</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.3 &times; 10<sup class='gt_super'>&minus;319</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #F0F5EF; color: #000000;">0.0060</td>
    ## <td class="gt_row gt_right">0.1396</td>
    ## <td class="gt_row gt_right">0.1396</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">1.653996e+02</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.2 &times; 10<sup class='gt_super'>&minus;37</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #DEF1D9; color: #000000;">0.0215</td>
    ## <td class="gt_row gt_right">1.7801</td>
    ## <td class="gt_row gt_right">1.7801</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">2.108706e+03</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F6F7; color: #000000;">&minus;0.0007</td>
    ## <td class="gt_row gt_right">0.0018</td>
    ## <td class="gt_row gt_right">0.0018</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">2.107968e+00</td>
    ## <td class="gt_row gt_right">0.15</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">H</td>
    ## <td class="gt_row gt_right" style="background-color: #F1F6EF; color: #000000;">0.0056</td>
    ## <td class="gt_row gt_right">0.1183</td>
    ## <td class="gt_row gt_right">0.1183</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">1.401947e+02</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.3 &times; 10<sup class='gt_super'>&minus;32</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">0.0000</td>
    ## <td class="gt_row gt_right">0.0000</td>
    ## <td class="gt_row gt_right">0.0000</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">5.586296e-03</td>
    ## <td class="gt_row gt_right">0.94</td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="8"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'pbkrtest.limit = 15360' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 15360)' or larger];
    ## but be warned that this may result in large computation time and memory use.
    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 15360' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 15360)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## <div id="tpzhtafrqn" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #tpzhtafrqn .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #tpzhtafrqn .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #tpzhtafrqn .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #tpzhtafrqn .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #tpzhtafrqn .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Contrast</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Estimate</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">SE</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Z ratio</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">p</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - grid<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #F4F0F4; color: #000000;">&minus;0.0021</td>
    ## <td class="gt_row gt_right">0.0007</td>
    ## <td class="gt_row gt_right">&minus;3.1344</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.009333831<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">envgeo - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F4F6; color: #000000;">&minus;0.0010</td>
    ## <td class="gt_row gt_right">0.0007</td>
    ## <td class="gt_row gt_right">&minus;1.5145</td>
    ## <td class="gt_row gt_right">0.428730954</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - trans<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #AF8DC3; color: #000000;">&minus;0.0223</td>
    ## <td class="gt_row gt_right">0.0007</td>
    ## <td class="gt_row gt_right">&minus;33.5784</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">grid - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #F4F6F4; color: #000000;">0.0011</td>
    ## <td class="gt_row gt_right">0.0007</td>
    ## <td class="gt_row gt_right">1.6199</td>
    ## <td class="gt_row gt_right">0.367294910</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">grid - trans<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #B99ACA; color: #000000;">&minus;0.0202</td>
    ## <td class="gt_row gt_right">0.0007</td>
    ## <td class="gt_row gt_right">&minus;30.4440</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">rand - trans<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #B493C6; color: #000000;">&minus;0.0213</td>
    ## <td class="gt_row gt_right">0.0007</td>
    ## <td class="gt_row gt_right">&minus;32.0639</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.01</td>
    ##     </tr>
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">**</sup> p &lt; 0.001</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

``` r
run_lmer(mmrr_ind, "env_ae", filepath = here(p4path, "MMRR_individual_enverr.csv"))
```

    ## <div id="fixjlqekcp" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #fixjlqekcp .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #fixjlqekcp .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #fixjlqekcp .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #fixjlqekcp .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #fixjlqekcp .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #fixjlqekcp .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #fixjlqekcp .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #fixjlqekcp .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #fixjlqekcp .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Predictors</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Fixed Effects</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Sum Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Mean Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">NumDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">DenDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">F value</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Pr(&gt;F)</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">nsamp</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;0.0001</td>
    ## <td class="gt_row gt_right">0.7452</td>
    ## <td class="gt_row gt_right">0.7452</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">1940.537011</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.0350</td>
    ## <td class="gt_row gt_right">0.0609</td>
    ## <td class="gt_row gt_right">0.0203</td>
    ## <td class="gt_row gt_right">3</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">52.882904</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">5.4 &times; 10<sup class='gt_super'>&minus;34</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #F4F6F4; color: #000000;">0.0017</td>
    ## <td class="gt_row gt_right">0.0105</td>
    ## <td class="gt_row gt_right">0.0105</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">27.383651</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.7 &times; 10<sup class='gt_super'>&minus;7</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #E7F3E3; color: #000000;">0.0095</td>
    ## <td class="gt_row gt_right">0.3480</td>
    ## <td class="gt_row gt_right">0.3480</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">906.190605</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.8 &times; 10<sup class='gt_super'>&minus;193</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F6F7; color: #000000;">&minus;0.0003</td>
    ## <td class="gt_row gt_right">0.0004</td>
    ## <td class="gt_row gt_right">0.0004</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">1.026801</td>
    ## <td class="gt_row gt_right">0.31</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">H</td>
    ## <td class="gt_row gt_right" style="background-color: #F5F6F4; color: #000000;">0.0014</td>
    ## <td class="gt_row gt_right">0.0078</td>
    ## <td class="gt_row gt_right">0.0078</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">20.317704</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">6.6 &times; 10<sup class='gt_super'>&minus;6</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F5F2F5; color: #000000;">&minus;0.0026</td>
    ## <td class="gt_row gt_right">0.0262</td>
    ## <td class="gt_row gt_right">0.0262</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3480K</td>
    ## <td class="gt_row gt_right">68.233041</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.6 &times; 10<sup class='gt_super'>&minus;16</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="8"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'pbkrtest.limit = 15360' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 15360)' or larger];
    ## but be warned that this may result in large computation time and memory use.
    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 15360' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 15360)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## <div id="tpzhtafrqn" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #tpzhtafrqn .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #tpzhtafrqn .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #tpzhtafrqn .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #tpzhtafrqn .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #tpzhtafrqn .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Contrast</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Estimate</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">SE</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Z ratio</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">p</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left">envgeo - grid</td>
    ## <td class="gt_row gt_right" style="background-color: #F4EFF4; color: #000000;">&minus;0.0005</td>
    ## <td class="gt_row gt_right">0.0004</td>
    ## <td class="gt_row gt_right">&minus;1.2011</td>
    ## <td class="gt_row gt_right">0.6261246</td></tr>
    ##     <tr><td class="gt_row gt_left">envgeo - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #F1F6F0; color: #000000;">0.0005</td>
    ## <td class="gt_row gt_right">0.0004</td>
    ## <td class="gt_row gt_right">1.0594</td>
    ## <td class="gt_row gt_right">0.7142981</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #BA9ACA; color: #000000;">&minus;0.0045</td>
    ## <td class="gt_row gt_right">0.0004</td>
    ## <td class="gt_row gt_right">&minus;10.1643</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.7 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">grid - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #EBF4E8; color: #000000;">0.0010</td>
    ## <td class="gt_row gt_right">0.0004</td>
    ## <td class="gt_row gt_right">2.2605</td>
    ## <td class="gt_row gt_right">0.1074219</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">grid - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #C6A9D2; color: #000000;">&minus;0.0040</td>
    ## <td class="gt_row gt_right">0.0004</td>
    ## <td class="gt_row gt_right">&minus;8.9632</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.7 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">rand - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #AF8DC3; color: #000000;">&minus;0.0050</td>
    ## <td class="gt_row gt_right">0.0004</td>
    ## <td class="gt_row gt_right">&minus;11.2237</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">8.0 &times; 10<sup class='gt_super'>&minus;15</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

### 1.1.3 Megaplots

``` r
MEGAPLOT(mmrr_ind, "ratio_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
MEGAPLOT(mmrr_ind, "ratio_err", colpal = "viridis", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
MEGAPLOT(mmrr_ind, "comboenv_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
MEGAPLOT(mmrr_ind, "geo_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->

## 1.2 Site sampling

### 1.2.1 Summary plots

``` r
mmrr_site <- format_mmrr(here(p3path, "mmrr_sitesampling_results.csv"))

# overall error
summary_hplot(mmrr_site, "ratio_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
summary_hplot(mmrr_site, "env_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
summary_hplot(mmrr_site, "geo_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
# bias
summary_hplot(mmrr_site, "ratio_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

``` r
summary_hplot(mmrr_site, "comboenv_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-5-5.png)<!-- -->

``` r
summary_hplot(mmrr_site, "geo_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-5-6.png)<!-- -->

### 1.2.2 Model summaries

``` r
run_lmer(mmrr_site, "ratio_ae", filepath = here(p4path, "MMRR_site_ratioerr.csv"))
```

    ## <div id="fixjlqekcp" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #fixjlqekcp .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #fixjlqekcp .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #fixjlqekcp .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #fixjlqekcp .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #fixjlqekcp .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #fixjlqekcp .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #fixjlqekcp .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #fixjlqekcp .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #fixjlqekcp .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Predictors</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Fixed Effects</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Sum Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Mean Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">NumDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">DenDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">F value</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Pr(&gt;F)</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">nsamp</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F5F6; color: #000000;">&minus;0.0054</td>
    ## <td class="gt_row gt_right">10.9296</td>
    ## <td class="gt_row gt_right">10.9296</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">774.65997483</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.4 &times; 10<sup class='gt_super'>&minus;163</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.2519</td>
    ## <td class="gt_row gt_right">1.7532</td>
    ## <td class="gt_row gt_right">0.8766</td>
    ## <td class="gt_row gt_right">2</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">62.13047517</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.6 &times; 10<sup class='gt_super'>&minus;27</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">K</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">0.0003</td>
    ## <td class="gt_row gt_right">0.0002</td>
    ## <td class="gt_row gt_right">0.0002</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">0.01177688</td>
    ## <td class="gt_row gt_right">0.910</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #F2F6F1; color: #000000;">0.0196</td>
    ## <td class="gt_row gt_right">0.8264</td>
    ## <td class="gt_row gt_right">0.8264</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">58.57599989</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.2 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F6F6; color: #000000;">&minus;0.0053</td>
    ## <td class="gt_row gt_right">0.0618</td>
    ## <td class="gt_row gt_right">0.0618</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">4.37735932</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.036<sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">H</td>
    ## <td class="gt_row gt_right" style="background-color: #F4F6F4; color: #000000;">0.0115</td>
    ## <td class="gt_row gt_right">0.2873</td>
    ## <td class="gt_row gt_right">0.2873</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">20.36485607</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">6.5 &times; 10<sup class='gt_super'>&minus;6</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F5F2F5; color: #000000;">&minus;0.0184</td>
    ## <td class="gt_row gt_right">0.7299</td>
    ## <td class="gt_row gt_right">0.7299</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">51.73067018</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">6.9 &times; 10<sup class='gt_super'>&minus;13</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="8"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##     <tr>
    ##       <td class="gt_footnote" colspan="8"><sup class="gt_footnote_marks gt_asterisk">**</sup> p &lt; 0.05</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'pbkrtest.limit = 8640' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 8640)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 8640' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 8640)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## <div id="tpzhtafrqn" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #tpzhtafrqn .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #tpzhtafrqn .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #tpzhtafrqn .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #tpzhtafrqn .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #tpzhtafrqn .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Contrast</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Estimate</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">SE</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Z ratio</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">p</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - equi<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #AF8DC3; color: #000000;">&minus;0.0340</td>
    ## <td class="gt_row gt_right">0.0031</td>
    ## <td class="gt_row gt_right">&minus;10.8514</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.7 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - rand<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #EEE2EE; color: #000000;">&minus;0.0101</td>
    ## <td class="gt_row gt_right">0.0031</td>
    ## <td class="gt_row gt_right">&minus;3.2163</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.7 &times; 10<sup class='gt_super'>&minus;3</sup><sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">equi - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #B5DCAF; color: #000000;">0.0239</td>
    ## <td class="gt_row gt_right">0.0031</td>
    ## <td class="gt_row gt_right">7.6351</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">8.6 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">**</sup> p &lt; 0.01</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

``` r
run_lmer(mmrr_site, "geo_ae", filepath = here(p4path, "MMRR_site_geoerr.csv"))
```

    ## <div id="fixjlqekcp" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #fixjlqekcp .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #fixjlqekcp .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #fixjlqekcp .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #fixjlqekcp .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #fixjlqekcp .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #fixjlqekcp .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #fixjlqekcp .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #fixjlqekcp .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #fixjlqekcp .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Predictors</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Fixed Effects</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Sum Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Mean Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">NumDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">DenDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">F value</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Pr(&gt;F)</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">nsamp</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;0.0011</td>
    ## <td class="gt_row gt_right">0.4225</td>
    ## <td class="gt_row gt_right">0.4225</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">131.605173</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.0 &times; 10<sup class='gt_super'>&minus;30</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #DDF1D8; color: #000000;">0.0702</td>
    ## <td class="gt_row gt_right">0.1566</td>
    ## <td class="gt_row gt_right">0.0783</td>
    ## <td class="gt_row gt_right">2</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">24.386603</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.7 &times; 10<sup class='gt_super'>&minus;11</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #EBF4E8; color: #000000;">0.0329</td>
    ## <td class="gt_row gt_right">2.3423</td>
    ## <td class="gt_row gt_right">2.3423</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">729.634765</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.5 &times; 10<sup class='gt_super'>&minus;154</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.1628</td>
    ## <td class="gt_row gt_right">57.2736</td>
    ## <td class="gt_row gt_right">57.2736</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">17841.045336</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F4F6; color: #000000;">&minus;0.0069</td>
    ## <td class="gt_row gt_right">0.1030</td>
    ## <td class="gt_row gt_right">0.1030</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">32.099252</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.5 &times; 10<sup class='gt_super'>&minus;8</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">H</td>
    ## <td class="gt_row gt_right" style="background-color: #F4F6F3; color: #000000;">0.0093</td>
    ## <td class="gt_row gt_right">0.1876</td>
    ## <td class="gt_row gt_right">0.1876</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">58.425672</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.3 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F7F5; color: #000000;">0.0036</td>
    ## <td class="gt_row gt_right">0.0275</td>
    ## <td class="gt_row gt_right">0.0275</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">8.564055</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.4 &times; 10<sup class='gt_super'>&minus;3</sup><sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="8"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##     <tr>
    ##       <td class="gt_footnote" colspan="8"><sup class="gt_footnote_marks gt_asterisk">**</sup> p &lt; 0.01</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'pbkrtest.limit = 8640' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 8640)' or larger];
    ## but be warned that this may result in large computation time and memory use.
    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 8640' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 8640)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## <div id="tpzhtafrqn" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #tpzhtafrqn .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #tpzhtafrqn .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #tpzhtafrqn .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #tpzhtafrqn .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #tpzhtafrqn .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Contrast</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Estimate</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">SE</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Z ratio</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">p</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - equi<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.0098</td>
    ## <td class="gt_row gt_right">0.0015</td>
    ## <td class="gt_row gt_right">6.5604</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.6 &times; 10<sup class='gt_super'>&minus;10</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">envgeo - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #ECF4EA; color: #000000;">0.0018</td>
    ## <td class="gt_row gt_right">0.0015</td>
    ## <td class="gt_row gt_right">1.2065</td>
    ## <td class="gt_row gt_right">0.4492773</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">equi - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #C4A7D1; color: #000000;">&minus;0.0080</td>
    ## <td class="gt_row gt_right">0.0015</td>
    ## <td class="gt_row gt_right">&minus;5.3539</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.6 &times; 10<sup class='gt_super'>&minus;7</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

``` r
run_lmer(mmrr_site, "env_ae", filepath = here(p4path, "MMRR_site_enverr.csv"))
```

    ## <div id="fixjlqekcp" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #fixjlqekcp .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #fixjlqekcp .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #fixjlqekcp .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #fixjlqekcp .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #fixjlqekcp .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #fixjlqekcp .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #fixjlqekcp .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #fixjlqekcp .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #fixjlqekcp .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Predictors</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Fixed Effects</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Sum Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Mean Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">NumDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">DenDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">F value</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Pr(&gt;F)</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">nsamp</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F5F6; color: #000000;">&minus;0.0021</td>
    ## <td class="gt_row gt_right">1.5946</td>
    ## <td class="gt_row gt_right">1.5946</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">907.15865165</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.4 &times; 10<sup class='gt_super'>&minus;189</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.0968</td>
    ## <td class="gt_row gt_right">0.2031</td>
    ## <td class="gt_row gt_right">0.1016</td>
    ## <td class="gt_row gt_right">2</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">57.77799916</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.2 &times; 10<sup class='gt_super'>&minus;25</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">K</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;0.0003</td>
    ## <td class="gt_row gt_right">0.0002</td>
    ## <td class="gt_row gt_right">0.0002</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">0.09364113</td>
    ## <td class="gt_row gt_right">0.760</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #F4F6F3; color: #000000;">0.0056</td>
    ## <td class="gt_row gt_right">0.0679</td>
    ## <td class="gt_row gt_right">0.0679</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">38.64283865</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">5.3 &times; 10<sup class='gt_super'>&minus;10</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F6F6; color: #000000;">&minus;0.0018</td>
    ## <td class="gt_row gt_right">0.0074</td>
    ## <td class="gt_row gt_right">0.0074</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">4.18983195</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.041<sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">H</td>
    ## <td class="gt_row gt_right" style="background-color: #F5F7F5; color: #000000;">0.0027</td>
    ## <td class="gt_row gt_right">0.0153</td>
    ## <td class="gt_row gt_right">0.0153</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">8.69527846</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.2 &times; 10<sup class='gt_super'>&minus;3</sup><sup class="gt_footnote_marks gt_asterisk">*</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F5F2F5; color: #000000;">&minus;0.0066</td>
    ## <td class="gt_row gt_right">0.0936</td>
    ## <td class="gt_row gt_right">0.0936</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.6290K</td>
    ## <td class="gt_row gt_right">53.23325588</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.2 &times; 10<sup class='gt_super'>&minus;13</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="8"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##     <tr>
    ##       <td class="gt_footnote" colspan="8"><sup class="gt_footnote_marks gt_asterisk">**</sup> p &lt; 0.05</td>
    ##     </tr>
    ##     <tr>
    ##       <td class="gt_footnote" colspan="8"><sup class="gt_footnote_marks gt_asterisk">*</sup> p &lt; 0.01</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'pbkrtest.limit = 8640' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 8640)' or larger];
    ## but be warned that this may result in large computation time and memory use.
    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 8640' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 8640)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## <div id="tpzhtafrqn" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #tpzhtafrqn .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #tpzhtafrqn .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #tpzhtafrqn .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #tpzhtafrqn .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #tpzhtafrqn .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Contrast</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Estimate</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">SE</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Z ratio</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">p</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - equi<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #AF8DC3; color: #000000;">&minus;0.0116</td>
    ## <td class="gt_row gt_right">0.0011</td>
    ## <td class="gt_row gt_right">&minus;10.4807</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.7 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - rand<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #EDE2EE; color: #000000;">&minus;0.0035</td>
    ## <td class="gt_row gt_right">0.0011</td>
    ## <td class="gt_row gt_right">&minus;3.1708</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">4.3 &times; 10<sup class='gt_super'>&minus;3</sup><sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">equi - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #B6DDB0; color: #000000;">0.0081</td>
    ## <td class="gt_row gt_right">0.0011</td>
    ## <td class="gt_row gt_right">7.3099</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">8.2 &times; 10<sup class='gt_super'>&minus;13</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">**</sup> p &lt; 0.01</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

### 1.2.3 Megaplots

``` r
MEGAPLOT(mmrr_site, "ratio_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
MEGAPLOT(mmrr_site, "ratio_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
MEGAPLOT(mmrr_site, "comboenv_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->

``` r
MEGAPLOT(mmrr_site, "geo_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-7-4.png)<!-- -->

# 2. GDM

## 2.1 Individual sampling

### 2.1.1 Summary plots

``` r
gdm_ind <- format_gdm(here(p3path, "gdm_indsampling_results.csv"))

# overall error
summary_hplot(gdm_ind, "ratio_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
summary_hplot(gdm_ind, "env_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
summary_hplot(gdm_ind, "geo_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->

``` r
# bias
summary_hplot(gdm_ind, "ratio_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->

``` r
summary_hplot(gdm_ind, "comboenv_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-8-5.png)<!-- -->

``` r
summary_hplot(gdm_ind, "geo_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-8-6.png)<!-- -->

### 2.1.2 Model summaries

``` r
run_lmer(gdm_ind, "ratio_ae", filepath = here(p4path, "GDM_individual_ratioerr.csv"))
```

    ## <div id="fixjlqekcp" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #fixjlqekcp .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #fixjlqekcp .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #fixjlqekcp .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #fixjlqekcp .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #fixjlqekcp .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #fixjlqekcp .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #fixjlqekcp .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #fixjlqekcp .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #fixjlqekcp .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Predictors</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Fixed Effects</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Sum Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Mean Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">NumDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">DenDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">F value</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Pr(&gt;F)</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">nsamp</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;0.0003</td>
    ## <td class="gt_row gt_right">9.4335</td>
    ## <td class="gt_row gt_right">9.4335</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">8.826316e+02</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.2 &times; 10<sup class='gt_super'>&minus;188</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.1099</td>
    ## <td class="gt_row gt_right">0.7896</td>
    ## <td class="gt_row gt_right">0.2632</td>
    ## <td class="gt_row gt_right">3</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">2.462613e+01</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">6.9 &times; 10<sup class='gt_super'>&minus;16</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #EAF4E8; color: #000000;">0.0231</td>
    ## <td class="gt_row gt_right">2.0533</td>
    ## <td class="gt_row gt_right">2.0533</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">1.921097e+02</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.0 &times; 10<sup class='gt_super'>&minus;43</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #B1DAAA; color: #000000;">0.0799</td>
    ## <td class="gt_row gt_right">24.4973</td>
    ## <td class="gt_row gt_right">24.4973</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">2.292047e+03</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;0.0003</td>
    ## <td class="gt_row gt_right">0.0002</td>
    ## <td class="gt_row gt_right">0.0002</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">2.244788e-02</td>
    ## <td class="gt_row gt_right">0.88</td></tr>
    ##     <tr><td class="gt_row gt_left">H</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F7F6; color: #000000;">0.0020</td>
    ## <td class="gt_row gt_right">0.0153</td>
    ## <td class="gt_row gt_right">0.0153</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">1.435012e+00</td>
    ## <td class="gt_row gt_right">0.23</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F5F3F5; color: #000000;">&minus;0.0057</td>
    ## <td class="gt_row gt_right">0.1237</td>
    ## <td class="gt_row gt_right">0.1237</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">1.157134e+01</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">6.7 &times; 10<sup class='gt_super'>&minus;4</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="8"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'pbkrtest.limit = 15355' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 15355)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 15355' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 15355)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## <div id="tpzhtafrqn" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #tpzhtafrqn .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #tpzhtafrqn .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #tpzhtafrqn .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #tpzhtafrqn .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #tpzhtafrqn .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Contrast</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Estimate</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">SE</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Z ratio</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">p</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - grid<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #D2B9DA; color: #000000;">&minus;0.0135</td>
    ## <td class="gt_row gt_right">0.0024</td>
    ## <td class="gt_row gt_right">&minus;5.7327</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">5.9 &times; 10<sup class='gt_super'>&minus;8</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #E7D5E8; color: #000000;">&minus;0.0097</td>
    ## <td class="gt_row gt_right">0.0024</td>
    ## <td class="gt_row gt_right">&minus;4.0895</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.5 &times; 10<sup class='gt_super'>&minus;4</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #AF8DC3; color: #000000;">&minus;0.0198</td>
    ## <td class="gt_row gt_right">0.0024</td>
    ## <td class="gt_row gt_right">&minus;8.3745</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">4.1 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">grid - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #EBF4E9; color: #000000;">0.0039</td>
    ## <td class="gt_row gt_right">0.0024</td>
    ## <td class="gt_row gt_right">1.6429</td>
    ## <td class="gt_row gt_right">0.35443857</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">grid - trans<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #EDE1EE; color: #000000;">&minus;0.0062</td>
    ## <td class="gt_row gt_right">0.0024</td>
    ## <td class="gt_row gt_right">&minus;2.6417</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.04111055<sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">rand - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #E6D2E7; color: #000000;">&minus;0.0101</td>
    ## <td class="gt_row gt_right">0.0024</td>
    ## <td class="gt_row gt_right">&minus;4.2845</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.1 &times; 10<sup class='gt_super'>&minus;4</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">**</sup> p &lt; 0.05</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

``` r
run_lmer(gdm_ind, "geo_ae", filepath = here(p4path, "GDM_individual_geoerr.csv"))
```

    ## <div id="fixjlqekcp" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #fixjlqekcp .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #fixjlqekcp .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #fixjlqekcp .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #fixjlqekcp .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #fixjlqekcp .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #fixjlqekcp .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #fixjlqekcp .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #fixjlqekcp .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #fixjlqekcp .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Predictors</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Fixed Effects</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Sum Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Mean Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">NumDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">DenDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">F value</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Pr(&gt;F)</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">nsamp</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;0.0012</td>
    ## <td class="gt_row gt_right">106.3637</td>
    ## <td class="gt_row gt_right">106.3637</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">2702.1766883</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.6274</td>
    ## <td class="gt_row gt_right">17.5095</td>
    ## <td class="gt_row gt_right">5.8365</td>
    ## <td class="gt_row gt_right">3</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">148.2767520</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.0 &times; 10<sup class='gt_super'>&minus;94</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #F5F6F4; color: #000000;">0.0235</td>
    ## <td class="gt_row gt_right">2.1250</td>
    ## <td class="gt_row gt_right">2.1250</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">53.9860038</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.1 &times; 10<sup class='gt_super'>&minus;13</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #F3EDF3; color: #000000;">&minus;0.0852</td>
    ## <td class="gt_row gt_right">27.8979</td>
    ## <td class="gt_row gt_right">27.8979</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">708.7482021</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.1 &times; 10<sup class='gt_super'>&minus;152</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F7F6; color: #000000;">0.0073</td>
    ## <td class="gt_row gt_right">0.2059</td>
    ## <td class="gt_row gt_right">0.2059</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">5.2302675</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.022<sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">H</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;0.0033</td>
    ## <td class="gt_row gt_right">0.0413</td>
    ## <td class="gt_row gt_right">0.0413</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">1.0486234</td>
    ## <td class="gt_row gt_right">0.310</td></tr>
    ##     <tr><td class="gt_row gt_left">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">0.0020</td>
    ## <td class="gt_row gt_right">0.0148</td>
    ## <td class="gt_row gt_right">0.0148</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">0.3751083</td>
    ## <td class="gt_row gt_right">0.540</td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="8"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##     <tr>
    ##       <td class="gt_footnote" colspan="8"><sup class="gt_footnote_marks gt_asterisk">**</sup> p &lt; 0.05</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'pbkrtest.limit = 15355' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 15355)' or larger];
    ## but be warned that this may result in large computation time and memory use.
    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 15355' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 15355)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## <div id="tpzhtafrqn" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #tpzhtafrqn .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #tpzhtafrqn .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #tpzhtafrqn .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #tpzhtafrqn .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #tpzhtafrqn .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Contrast</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Estimate</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">SE</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Z ratio</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">p</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - grid<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.0899</td>
    ## <td class="gt_row gt_right">0.0045</td>
    ## <td class="gt_row gt_right">19.8550</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #E9F4E7; color: #000000;">0.0204</td>
    ## <td class="gt_row gt_right">0.0045</td>
    ## <td class="gt_row gt_right">4.5115</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.8 &times; 10<sup class='gt_super'>&minus;5</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #D0EBCA; color: #000000;">0.0493</td>
    ## <td class="gt_row gt_right">0.0045</td>
    ## <td class="gt_row gt_right">10.8910</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">4.0 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">grid - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #C8ADD4; color: #000000;">&minus;0.0695</td>
    ## <td class="gt_row gt_right">0.0045</td>
    ## <td class="gt_row gt_right">&minus;15.3432</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">grid - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #E9D7E9; color: #000000;">&minus;0.0406</td>
    ## <td class="gt_row gt_right">0.0045</td>
    ## <td class="gt_row gt_right">&minus;8.9659</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.7 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">rand - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #E4F3E0; color: #000000;">0.0289</td>
    ## <td class="gt_row gt_right">0.0045</td>
    ## <td class="gt_row gt_right">6.3789</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.1 &times; 10<sup class='gt_super'>&minus;9</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

``` r
run_lmer(gdm_ind, "env_ae", filepath = here(p4path, "GDM_individual_enverr.csv"))
```

    ## <div id="fixjlqekcp" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #fixjlqekcp .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #fixjlqekcp .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #fixjlqekcp .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #fixjlqekcp .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #fixjlqekcp .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #fixjlqekcp .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #fixjlqekcp .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #fixjlqekcp .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #fixjlqekcp .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Predictors</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Fixed Effects</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Sum Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Mean Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">NumDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">DenDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">F value</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Pr(&gt;F)</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">nsamp</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;0.0003</td>
    ## <td class="gt_row gt_right">5.3338</td>
    ## <td class="gt_row gt_right">5.3338</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">1716.7973944</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.0999</td>
    ## <td class="gt_row gt_right">0.2357</td>
    ## <td class="gt_row gt_right">0.0786</td>
    ## <td class="gt_row gt_right">3</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">25.2933088</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.6 &times; 10<sup class='gt_super'>&minus;16</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F5F6; color: #000000;">&minus;0.0034</td>
    ## <td class="gt_row gt_right">0.0432</td>
    ## <td class="gt_row gt_right">0.0432</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">13.8946308</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.9 &times; 10<sup class='gt_super'>&minus;4</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #F3F6F2; color: #000000;">0.0062</td>
    ## <td class="gt_row gt_right">0.1486</td>
    ## <td class="gt_row gt_right">0.1486</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">47.8286525</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">4.8 &times; 10<sup class='gt_super'>&minus;12</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #F4F6F4; color: #000000;">0.0043</td>
    ## <td class="gt_row gt_right">0.0716</td>
    ## <td class="gt_row gt_right">0.0716</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">23.0482092</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.6 &times; 10<sup class='gt_super'>&minus;6</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">H</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F7F6; color: #000000;">0.0009</td>
    ## <td class="gt_row gt_right">0.0030</td>
    ## <td class="gt_row gt_right">0.0030</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">0.9623282</td>
    ## <td class="gt_row gt_right">0.33</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F4F6; color: #000000;">&minus;0.0040</td>
    ## <td class="gt_row gt_right">0.0619</td>
    ## <td class="gt_row gt_right">0.0619</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">15.3430K</td>
    ## <td class="gt_row gt_right">19.9363228</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">8.1 &times; 10<sup class='gt_super'>&minus;6</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="8"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'pbkrtest.limit = 15355' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 15355)' or larger];
    ## but be warned that this may result in large computation time and memory use.
    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 15355' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 15355)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## <div id="tpzhtafrqn" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #tpzhtafrqn .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #tpzhtafrqn .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #tpzhtafrqn .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #tpzhtafrqn .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #tpzhtafrqn .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Contrast</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Estimate</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">SE</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Z ratio</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">p</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left">envgeo - grid</td>
    ## <td class="gt_row gt_right" style="background-color: #EFE6F0; color: #000000;">&minus;0.0026</td>
    ## <td class="gt_row gt_right">0.0013</td>
    ## <td class="gt_row gt_right">&minus;2.0550</td>
    ## <td class="gt_row gt_right">0.1681013</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #E9D8EA; color: #000000;">&minus;0.0046</td>
    ## <td class="gt_row gt_right">0.0013</td>
    ## <td class="gt_row gt_right">&minus;3.6538</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.5 &times; 10<sup class='gt_super'>&minus;3</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - trans<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #AF8DC3; color: #000000;">&minus;0.0106</td>
    ## <td class="gt_row gt_right">0.0013</td>
    ## <td class="gt_row gt_right">&minus;8.3556</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">4.3 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">grid - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #F1EAF1; color: #000000;">&minus;0.0020</td>
    ## <td class="gt_row gt_right">0.0013</td>
    ## <td class="gt_row gt_right">&minus;1.5990</td>
    ## <td class="gt_row gt_right">0.3791273</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">grid - trans<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #CBAFD5; color: #000000;">&minus;0.0080</td>
    ## <td class="gt_row gt_right">0.0013</td>
    ## <td class="gt_row gt_right">&minus;6.3009</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.8 &times; 10<sup class='gt_super'>&minus;9</sup><sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">rand - trans<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #E0CBE3; color: #000000;">&minus;0.0060</td>
    ## <td class="gt_row gt_right">0.0013</td>
    ## <td class="gt_row gt_right">&minus;4.7013</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.5 &times; 10<sup class='gt_super'>&minus;5</sup><sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.01</td>
    ##     </tr>
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">**</sup> p &lt; 0.001</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

### 2.1.3 Megaplots

``` r
MEGAPLOT(gdm_ind, "ratio_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
MEGAPLOT(gdm_ind, "ratio_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

``` r
MEGAPLOT(gdm_ind, "comboenv_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->

``` r
MEGAPLOT(gdm_ind, "geo_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-10-4.png)<!-- -->

### 2.1.4 Prop NA

Confirm that the distribution of NAs is as expected and the proportions
are small

``` r
MEGAPLOT(gdm_ind, "geo_coeff", aggfunc = "prop_na", colpal = "mako")
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

## 2.2 Site sampling

### 2.2.1 Summary plots

``` r
gdm_site <- format_gdm(here(p3path, "gdm_sitesampling_results.csv"))

# overall error
summary_hplot(gdm_site, "ratio_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
summary_hplot(gdm_site, "env_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

``` r
summary_hplot(gdm_site, "geo_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-12-3.png)<!-- -->

``` r
# bias
summary_hplot(gdm_site, "ratio_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-12-4.png)<!-- -->

``` r
summary_hplot(gdm_site, "comboenv_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-12-5.png)<!-- -->

``` r
summary_hplot(gdm_site, "geo_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-12-6.png)<!-- -->

### 2.2.2 Model summaries

``` r
run_lmer(gdm_site, "ratio_ae", filepath = here(p4path, "GDM_site_ratioerr.csv"))
```

    ## <div id="fixjlqekcp" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #fixjlqekcp .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #fixjlqekcp .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #fixjlqekcp .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #fixjlqekcp .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #fixjlqekcp .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #fixjlqekcp .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #fixjlqekcp .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #fixjlqekcp .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #fixjlqekcp .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Predictors</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Fixed Effects</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Sum Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Mean Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">NumDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">DenDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">F value</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Pr(&gt;F)</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">nsamp</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F6F6; color: #000000;">&minus;0.0064</td>
    ## <td class="gt_row gt_right">14.4681</td>
    ## <td class="gt_row gt_right">14.4681</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">5.407432e+02</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">5.9 &times; 10<sup class='gt_super'>&minus;116</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.3179</td>
    ## <td class="gt_row gt_right">8.3211</td>
    ## <td class="gt_row gt_right">4.1605</td>
    ## <td class="gt_row gt_right">2</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">1.554993e+02</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">5.0 &times; 10<sup class='gt_super'>&minus;67</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">K</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">0.0000</td>
    ## <td class="gt_row gt_right">0.0000</td>
    ## <td class="gt_row gt_right">0.0000</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">2.777680e-05</td>
    ## <td class="gt_row gt_right">1.00</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #F3F6F2; color: #000000;">0.0231</td>
    ## <td class="gt_row gt_right">1.1117</td>
    ## <td class="gt_row gt_right">1.1117</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">4.155104e+01</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.2 &times; 10<sup class='gt_super'>&minus;10</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F6F7; color: #000000;">&minus;0.0049</td>
    ## <td class="gt_row gt_right">0.0503</td>
    ## <td class="gt_row gt_right">0.0503</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">1.879243e+00</td>
    ## <td class="gt_row gt_right">0.17</td></tr>
    ##     <tr><td class="gt_row gt_left">H</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">0.0007</td>
    ## <td class="gt_row gt_right">0.0009</td>
    ## <td class="gt_row gt_right">0.0009</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">3.448557e-02</td>
    ## <td class="gt_row gt_right">0.85</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F4F6; color: #000000;">&minus;0.0137</td>
    ## <td class="gt_row gt_right">0.3905</td>
    ## <td class="gt_row gt_right">0.3905</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">1.459672e+01</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.3 &times; 10<sup class='gt_super'>&minus;4</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="8"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'pbkrtest.limit = 8361' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 8361)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 8361' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 8361)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## <div id="tpzhtafrqn" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #tpzhtafrqn .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #tpzhtafrqn .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #tpzhtafrqn .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #tpzhtafrqn .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #tpzhtafrqn .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Contrast</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Estimate</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">SE</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Z ratio</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">p</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - equi<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #AF8DC3; color: #000000;">&minus;0.0760</td>
    ## <td class="gt_row gt_right">0.0044</td>
    ## <td class="gt_row gt_right">&minus;17.3557</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #ECDEEC; color: #000000;">&minus;0.0272</td>
    ## <td class="gt_row gt_right">0.0044</td>
    ## <td class="gt_row gt_right">&minus;6.1624</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.1 &times; 10<sup class='gt_super'>&minus;9</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">equi - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #C0E2BA; color: #000000;">0.0488</td>
    ## <td class="gt_row gt_right">0.0044</td>
    ## <td class="gt_row gt_right">11.1986</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">7.3 &times; 10<sup class='gt_super'>&minus;15</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

``` r
run_lmer(gdm_site, "geo_ae", filepath = here(p4path, "GDM_site_geoerr.csv"))
```

    ## <div id="fixjlqekcp" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #fixjlqekcp .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #fixjlqekcp .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #fixjlqekcp .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #fixjlqekcp .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #fixjlqekcp .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #fixjlqekcp .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #fixjlqekcp .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #fixjlqekcp .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #fixjlqekcp .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Predictors</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Fixed Effects</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Sum Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Mean Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">NumDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">DenDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">F value</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Pr(&gt;F)</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">nsamp</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F6F7; color: #000000;">&minus;0.0102</td>
    ## <td class="gt_row gt_right">37.3428</td>
    ## <td class="gt_row gt_right">37.3428</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3501K</td>
    ## <td class="gt_row gt_right">147.9671384</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">9.3 &times; 10<sup class='gt_super'>&minus;34</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">1.1474</td>
    ## <td class="gt_row gt_right">145.0126</td>
    ## <td class="gt_row gt_right">72.5063</td>
    ## <td class="gt_row gt_right">2</td>
    ## <td class="gt_row gt_right">8.3501K</td>
    ## <td class="gt_row gt_right">287.2989940</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.2 &times; 10<sup class='gt_super'>&minus;121</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #EFF5EE; color: #000000;">0.1451</td>
    ## <td class="gt_row gt_right">43.9976</td>
    ## <td class="gt_row gt_right">43.9976</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">174.3360187</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.1 &times; 10<sup class='gt_super'>&minus;39</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #D1ECCB; color: #000000;">0.6239</td>
    ## <td class="gt_row gt_right">813.6420</td>
    ## <td class="gt_row gt_right">813.6420</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">3223.9742765</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F6; color: #000000;">0.0080</td>
    ## <td class="gt_row gt_right">0.1332</td>
    ## <td class="gt_row gt_right">0.1332</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">0.5279855</td>
    ## <td class="gt_row gt_right">0.470</td></tr>
    ##     <tr><td class="gt_row gt_left">H</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">0.0041</td>
    ## <td class="gt_row gt_right">0.0344</td>
    ## <td class="gt_row gt_right">0.0344</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">0.1361223</td>
    ## <td class="gt_row gt_right">0.710</td></tr>
    ##     <tr><td class="gt_row gt_left">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F7F6; color: #000000;">0.0201</td>
    ## <td class="gt_row gt_right">0.8420</td>
    ## <td class="gt_row gt_right">0.8420</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3501K</td>
    ## <td class="gt_row gt_right">3.3362023</td>
    ## <td class="gt_row gt_right">0.068</td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="8"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'pbkrtest.limit = 8361' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 8361)' or larger];
    ## but be warned that this may result in large computation time and memory use.
    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 8361' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 8361)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## <div id="tpzhtafrqn" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #tpzhtafrqn .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #tpzhtafrqn .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #tpzhtafrqn .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #tpzhtafrqn .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #tpzhtafrqn .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Contrast</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Estimate</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">SE</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Z ratio</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">p</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - equi<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.2882</td>
    ## <td class="gt_row gt_right">0.0134</td>
    ## <td class="gt_row gt_right">21.4290</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">envgeo - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #F2F6F1; color: #000000;">0.0229</td>
    ## <td class="gt_row gt_right">0.0136</td>
    ## <td class="gt_row gt_right">1.6897</td>
    ## <td class="gt_row gt_right">0.2091145</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">equi - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #B898C9; color: #000000;">&minus;0.2652</td>
    ## <td class="gt_row gt_right">0.0134</td>
    ## <td class="gt_row gt_right">&minus;19.8333</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

``` r
run_lmer(gdm_site, "env_ae", filepath = here(p4path, "GDM_site_enverr.csv"))
```

    ## <div id="fixjlqekcp" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #fixjlqekcp .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #fixjlqekcp .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #fixjlqekcp .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #fixjlqekcp .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #fixjlqekcp .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #fixjlqekcp .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #fixjlqekcp .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #fixjlqekcp .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #fixjlqekcp .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #fixjlqekcp .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #fixjlqekcp .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #fixjlqekcp .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #fixjlqekcp .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #fixjlqekcp .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #fixjlqekcp .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #fixjlqekcp .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Predictors</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Fixed Effects</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Sum Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Mean Sq</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">NumDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">DenDF</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">F value</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Pr(&gt;F)</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">nsamp</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F6F6; color: #000000;">&minus;0.0054</td>
    ## <td class="gt_row gt_right">10.2393</td>
    ## <td class="gt_row gt_right">10.2393</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">625.37735767</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.7 &times; 10<sup class='gt_super'>&minus;133</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.2722</td>
    ## <td class="gt_row gt_right">1.5358</td>
    ## <td class="gt_row gt_right">0.7679</td>
    ## <td class="gt_row gt_right">2</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">46.90111271</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">5.6 &times; 10<sup class='gt_super'>&minus;21</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">K</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">0.0003</td>
    ## <td class="gt_row gt_right">0.0002</td>
    ## <td class="gt_row gt_right">0.0002</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">0.01089743</td>
    ## <td class="gt_row gt_right">0.920</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #F2F6F1; color: #000000;">0.0236</td>
    ## <td class="gt_row gt_right">1.1665</td>
    ## <td class="gt_row gt_right">1.1665</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">71.24432149</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.7 &times; 10<sup class='gt_super'>&minus;17</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">0.0010</td>
    ## <td class="gt_row gt_right">0.0022</td>
    ## <td class="gt_row gt_right">0.0022</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">0.13328216</td>
    ## <td class="gt_row gt_right">0.720</td></tr>
    ##     <tr><td class="gt_row gt_left">H</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F7F6; color: #000000;">0.0051</td>
    ## <td class="gt_row gt_right">0.0544</td>
    ## <td class="gt_row gt_right">0.0544</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">3.32312411</td>
    ## <td class="gt_row gt_right">0.068</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F5F3F5; color: #000000;">&minus;0.0149</td>
    ## <td class="gt_row gt_right">0.4610</td>
    ## <td class="gt_row gt_right">0.4610</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">8.3500K</td>
    ## <td class="gt_row gt_right">28.15629246</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.1 &times; 10<sup class='gt_super'>&minus;7</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="8"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'pbkrtest.limit = 8361' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 8361)' or larger];
    ## but be warned that this may result in large computation time and memory use.
    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 8361' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 8361)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## <div id="tpzhtafrqn" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
    ##   <style>html {
    ##   font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
    ## }
    ## 
    ## #tpzhtafrqn .gt_table {
    ##   display: table;
    ##   border-collapse: collapse;
    ##   margin-left: auto;
    ##   margin-right: auto;
    ##   color: #333333;
    ##   font-size: 16px;
    ##   font-weight: normal;
    ##   font-style: normal;
    ##   background-color: #FFFFFF;
    ##   width: auto;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #A8A8A8;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #A8A8A8;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_heading {
    ##   background-color: #FFFFFF;
    ##   text-align: center;
    ##   border-bottom-color: #FFFFFF;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_title {
    ##   color: #333333;
    ##   font-size: 125%;
    ##   font-weight: initial;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-color: #FFFFFF;
    ##   border-bottom-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_subtitle {
    ##   color: #333333;
    ##   font-size: 85%;
    ##   font-weight: initial;
    ##   padding-top: 0;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-color: #FFFFFF;
    ##   border-top-width: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_bottom_border {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_headings {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_col_heading {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 6px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: normal;
    ##   text-transform: inherit;
    ##   padding-top: 0;
    ##   padding-bottom: 0;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:first-child {
    ##   padding-left: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner_outer:last-child {
    ##   padding-right: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_column_spanner {
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: bottom;
    ##   padding-top: 5px;
    ##   padding-bottom: 5px;
    ##   overflow-x: hidden;
    ##   display: inline-block;
    ##   width: 100%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_group_heading {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_empty_group_heading {
    ##   padding: 0.5px;
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   vertical-align: middle;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :first-child {
    ##   margin-top: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_from_md > :last-child {
    ##   margin-bottom: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   margin: 10px;
    ##   border-top-style: solid;
    ##   border-top-width: 1px;
    ##   border-top-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 1px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 1px;
    ##   border-right-color: #D3D3D3;
    ##   vertical-align: middle;
    ##   overflow-x: hidden;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_stub_row_group {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   font-size: 100%;
    ##   font-weight: initial;
    ##   text-transform: inherit;
    ##   border-right-style: solid;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   vertical-align: top;
    ## }
    ## 
    ## #tpzhtafrqn .gt_row_group_first td {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row {
    ##   border-top-style: solid;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_summary_row.thick {
    ##   border-top-width: 2px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_last_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_grand_summary_row {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   text-transform: inherit;
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_first_grand_summary_row {
    ##   padding-top: 8px;
    ##   padding-bottom: 8px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ##   border-top-style: double;
    ##   border-top-width: 6px;
    ##   border-top-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_striped {
    ##   background-color: rgba(128, 128, 128, 0.05);
    ## }
    ## 
    ## #tpzhtafrqn .gt_table_body {
    ##   border-top-style: solid;
    ##   border-top-width: 2px;
    ##   border-top-color: #D3D3D3;
    ##   border-bottom-style: solid;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote {
    ##   margin: 0px;
    ##   font-size: 90%;
    ##   padding-left: 4px;
    ##   padding-right: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenotes {
    ##   color: #333333;
    ##   background-color: #FFFFFF;
    ##   border-bottom-style: none;
    ##   border-bottom-width: 2px;
    ##   border-bottom-color: #D3D3D3;
    ##   border-left-style: none;
    ##   border-left-width: 2px;
    ##   border-left-color: #D3D3D3;
    ##   border-right-style: none;
    ##   border-right-width: 2px;
    ##   border-right-color: #D3D3D3;
    ## }
    ## 
    ## #tpzhtafrqn .gt_sourcenote {
    ##   font-size: 90%;
    ##   padding-top: 4px;
    ##   padding-bottom: 4px;
    ##   padding-left: 5px;
    ##   padding-right: 5px;
    ## }
    ## 
    ## #tpzhtafrqn .gt_left {
    ##   text-align: left;
    ## }
    ## 
    ## #tpzhtafrqn .gt_center {
    ##   text-align: center;
    ## }
    ## 
    ## #tpzhtafrqn .gt_right {
    ##   text-align: right;
    ##   font-variant-numeric: tabular-nums;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_normal {
    ##   font-weight: normal;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_bold {
    ##   font-weight: bold;
    ## }
    ## 
    ## #tpzhtafrqn .gt_font_italic {
    ##   font-style: italic;
    ## }
    ## 
    ## #tpzhtafrqn .gt_super {
    ##   font-size: 65%;
    ## }
    ## 
    ## #tpzhtafrqn .gt_two_val_uncert {
    ##   display: inline-block;
    ##   line-height: 1em;
    ##   text-align: right;
    ##   font-size: 60%;
    ##   vertical-align: -0.25em;
    ##   margin-left: 0.1em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_footnote_marks {
    ##   font-style: italic;
    ##   font-weight: normal;
    ##   font-size: 75%;
    ##   vertical-align: 0.4em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_asterisk {
    ##   font-size: 100%;
    ##   vertical-align: 0;
    ## }
    ## 
    ## #tpzhtafrqn .gt_slash_mark {
    ##   font-size: 0.7em;
    ##   line-height: 0.7em;
    ##   vertical-align: 0.15em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_numerator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: 0.45em;
    ## }
    ## 
    ## #tpzhtafrqn .gt_fraction_denominator {
    ##   font-size: 0.6em;
    ##   line-height: 0.6em;
    ##   vertical-align: -0.05em;
    ## }
    ## </style>
    ##   <table class="gt_table">
    ##   
    ##   <thead class="gt_col_headings">
    ##     <tr>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">Contrast</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Estimate</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">SE</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">Z ratio</th>
    ##       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">p</th>
    ##     </tr>
    ##   </thead>
    ##   <tbody class="gt_table_body">
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - equi<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #AF8DC3; color: #000000;">&minus;0.0307</td>
    ## <td class="gt_row gt_right">0.0034</td>
    ## <td class="gt_row gt_right">&minus;8.9621</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.1 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #BE9FCD; color: #000000;">&minus;0.0267</td>
    ## <td class="gt_row gt_right">0.0035</td>
    ## <td class="gt_row gt_right">&minus;7.7294</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">5.0 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">equi - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #EFF5EE; color: #000000;">0.0040</td>
    ## <td class="gt_row gt_right">0.0034</td>
    ## <td class="gt_row gt_right">1.1684</td>
    ## <td class="gt_row gt_right">0.4721507</td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.001</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

### 2.2.3 Megaplots

``` r
MEGAPLOT(gdm_site, "ratio_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
MEGAPLOT(gdm_site, "ratio_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

``` r
MEGAPLOT(gdm_site, "ratio_ae", colpal = "viridis", direction = -1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-14-3.png)<!-- -->

``` r
MEGAPLOT(gdm_site, "ratio_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-14-4.png)<!-- -->

``` r
MEGAPLOT(gdm_site, "comboenv_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-14-5.png)<!-- -->

``` r
MEGAPLOT(gdm_site, "geo_err", divergent = TRUE)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-14-6.png)<!-- -->

### 2.2.4 Prop NA

Confirm that the distribution of NAs is as expected and the proportions
are small

``` r
MEGAPLOT(gdm_site, "geo_coeff", aggfunc = "prop_na", colpal = "mako")
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

# 3. Comparison of error distirbutions for MMRR and GDM

``` r
mmrr_ind %>%
  filter(sampstrat != "full") %>%
  ggplot(aes(x = ratio_ae, fill = m, colour = m)) +
  geom_density(alpha = 0.5) +
  theme_few() +
  scale_fill_viridis_d(direction = -1, end = 0.7, begin = 0.3, option = "mako") +
  scale_color_viridis_d(direction = -1, end = 0.7, begin = 0.3, option = "mako")
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
gdm_ind %>%
  filter(sampstrat != "full") %>%
  ggplot(aes(x = ratio_ae, fill = m, colour = m)) +
  geom_density(alpha = 0.5) +
  theme_few() +
  scale_fill_viridis_d(direction = -1, end = 0.7, begin = 0.3, option = "mako") +
  scale_color_viridis_d(direction = -1, end = 0.7, begin = 0.3, option = "mako")
```

    ## Warning: Removed 5 rows containing non-finite values (stat_density).

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

# 4. Compare results of MMRR and GDM

``` r
plot(mmrr_ind$geo_coeff, gdm_ind$geo_coeff, col = gdm_ind$m)
legend("topleft", c("0.25", "1"), col = c("black", "red"), pch = 1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
plot(mmrr_ind$comboenv_coeff, gdm_ind$comboenv_coeff)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->

``` r
df <- data.frame(mmrr_ind, 
                 geo_mg = gdm_ind$geo_coeff/mmrr_ind$geo_coeff, 
                 env_mg = gdm_ind$comboenv_coeff/mmrr_ind$comboenv_coeff,
                 ratio_mg = gdm_ind$ratio/mmrr_ind$ratio )
```

``` r
summary_hplot(df, "geo_mg")
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
summary_hplot(df, "env_mg")
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

``` r
summary_hplot(df, "ratio_mg")
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-18-3.png)<!-- -->

``` r
MEGAPLOT(df, "geo_mg")
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
MEGAPLOT(df, "env_mg")
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

``` r
MEGAPLOT(df, "ratio_mg")
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-19-3.png)<!-- -->

``` r
plot(mmrr_site$geo_coeff, gdm_site$geo_coeff, col = gdm_site$m)
legend("topleft", c("0.25", "1"), col = c("black", "red"), pch = 1)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
plot(mmrr_site$comboenv_coeff, gdm_site$comboenv_coeff)
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->

``` r
df <- data.frame(mmrr_site, 
                 geo_mg = gdm_site$geo_coeff/mmrr_site$geo_coeff, 
                 env_mg = gdm_site$comboenv_coeff/mmrr_site$comboenv_coeff,
                 ratio_mg = gdm_site$ratio/mmrr_site$ratio )
```

``` r
summary_hplot(df, "geo_mg")
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
summary_hplot(df, "env_mg")
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-21-2.png)<!-- -->

``` r
summary_hplot(df, "ratio_mg")
```

![](analysis_IBDIBE_files/figure-gfm/unnamed-chunk-21-3.png)<!-- -->
