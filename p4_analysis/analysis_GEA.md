GEA Analysis
================

# 1. LFMM

## 1.1 Individual sampling

### 1.1.1 Summary plots

### K

``` r
walk(K_plots, grid.arrange)
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-3-4.png)<!-- -->

### TPRCOMBO

``` r
walk(TPR_plots, grid.arrange)
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->

### FDRCOMBO

``` r
walk(FDR_plots, grid.arrange)
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->
\### TOTALN

``` r
walk(TOTALN_plots, grid.arrange)
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->

### 1.1.2 Model summaries

``` r
lfmm_ind <- 
  lfmm_ind %>% 
  filter(method == "ridge") %>%
  filter(K_selection == "tracy.widom")
```

``` r
run_lmer(lfmm_ind, "TPRCOMBO", filepath = here(p4path, "LFMM_individual_TPR.csv"))
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
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">0.0003</td>
    ## <td class="gt_row gt_right">53.9590</td>
    ## <td class="gt_row gt_right">53.9590</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8680K</td>
    ## <td class="gt_row gt_right">1591.4994</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #CCE9C6; color: #000000;">0.1034</td>
    ## <td class="gt_row gt_right">20.1079</td>
    ## <td class="gt_row gt_right">6.7026</td>
    ## <td class="gt_row gt_right">3</td>
    ## <td class="gt_row gt_right">122.8680K</td>
    ## <td class="gt_row gt_right">197.6912</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">6.5 &times; 10<sup class='gt_super'>&minus;128</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #E9F4E5; color: #000000;">0.0437</td>
    ## <td class="gt_row gt_right">58.7945</td>
    ## <td class="gt_row gt_right">58.7945</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8680K</td>
    ## <td class="gt_row gt_right">1734.1216</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #F1F6F0; color: #000000;">0.0182</td>
    ## <td class="gt_row gt_right">10.2243</td>
    ## <td class="gt_row gt_right">10.2243</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8680K</td>
    ## <td class="gt_row gt_right">301.5614</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.8 &times; 10<sup class='gt_super'>&minus;67</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #B4DCAE; color: #000000;">0.1277</td>
    ## <td class="gt_row gt_right">501.1818</td>
    ## <td class="gt_row gt_right">501.1818</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8680K</td>
    ## <td class="gt_row gt_right">14782.1597</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">H</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.1805</td>
    ## <td class="gt_row gt_right">1.0004K</td>
    ## <td class="gt_row gt_right">1.0004K</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8680K</td>
    ## <td class="gt_row gt_right">29505.2496</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F2ECF2; color: #000000;">&minus;0.0276</td>
    ## <td class="gt_row gt_right">23.4670</td>
    ## <td class="gt_row gt_right">23.4670</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8680K</td>
    ## <td class="gt_row gt_right">692.1508</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">4.0 &times; 10<sup class='gt_super'>&minus;152</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
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
    ## To enable adjustments, add the argument 'pbkrtest.limit = 122880' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 122880)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 122880' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 122880)' or larger];
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
    ## <td class="gt_row gt_right" style="background-color: #F5F2F5; color: #000000;">&minus;0.0021</td>
    ## <td class="gt_row gt_right">0.0015</td>
    ## <td class="gt_row gt_right">&minus;1.4324</td>
    ## <td class="gt_row gt_right">0.47901837</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #F0F5EE; color: #000000;">0.0038</td>
    ## <td class="gt_row gt_right">0.0015</td>
    ## <td class="gt_row gt_right">2.5799</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.04857191<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - trans<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #8BC686; color: #000000;">0.0297</td>
    ## <td class="gt_row gt_right">0.0015</td>
    ## <td class="gt_row gt_right">19.9876</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">grid - rand<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #ECF4E9; color: #000000;">0.0060</td>
    ## <td class="gt_row gt_right">0.0015</td>
    ## <td class="gt_row gt_right">4.0123</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.5 &times; 10<sup class='gt_super'>&minus;4</sup><sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">grid - trans<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.0318</td>
    ## <td class="gt_row gt_right">0.0015</td>
    ## <td class="gt_row gt_right">21.4199</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">rand - trans<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #A1D19B; color: #000000;">0.0259</td>
    ## <td class="gt_row gt_right">0.0015</td>
    ## <td class="gt_row gt_right">17.4076</td>
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
run_lmer(lfmm_ind, "FDRCOMBO", filepath = here(p4path, "LFMM_individual_FDR.csv"))
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
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;0.0020</td>
    ## <td class="gt_row gt_right">2.3956K</td>
    ## <td class="gt_row gt_right">2.3956K</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8680K</td>
    ## <td class="gt_row gt_right">11734.92527</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.9315</td>
    ## <td class="gt_row gt_right">75.2508</td>
    ## <td class="gt_row gt_right">25.0836</td>
    ## <td class="gt_row gt_right">3</td>
    ## <td class="gt_row gt_right">122.8680K</td>
    ## <td class="gt_row gt_right">122.87478</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.8 &times; 10<sup class='gt_super'>&minus;79</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F4F6; color: #000000;">&minus;0.0363</td>
    ## <td class="gt_row gt_right">40.4242</td>
    ## <td class="gt_row gt_right">40.4242</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8680K</td>
    ## <td class="gt_row gt_right">198.02232</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">6.1 &times; 10<sup class='gt_super'>&minus;45</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #F5F3F5; color: #000000;">&minus;0.0518</td>
    ## <td class="gt_row gt_right">82.3529</td>
    ## <td class="gt_row gt_right">82.3529</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8680K</td>
    ## <td class="gt_row gt_right">403.41497</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.4 &times; 10<sup class='gt_super'>&minus;89</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #F5F3F5; color: #000000;">&minus;0.0486</td>
    ## <td class="gt_row gt_right">72.5919</td>
    ## <td class="gt_row gt_right">72.5919</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8680K</td>
    ## <td class="gt_row gt_right">355.59931</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.3 &times; 10<sup class='gt_super'>&minus;79</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">H</td>
    ## <td class="gt_row gt_right" style="background-color: #F5F7F5; color: #000000;">0.0260</td>
    ## <td class="gt_row gt_right">20.8290</td>
    ## <td class="gt_row gt_right">20.8290</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8680K</td>
    ## <td class="gt_row gt_right">102.03316</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">5.6 &times; 10<sup class='gt_super'>&minus;24</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;0.0048</td>
    ## <td class="gt_row gt_right">0.6950</td>
    ## <td class="gt_row gt_right">0.6950</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8680K</td>
    ## <td class="gt_row gt_right">3.40472</td>
    ## <td class="gt_row gt_right">0.065</td></tr>
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
    ## To enable adjustments, add the argument 'pbkrtest.limit = 122880' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 122880)' or larger];
    ## but be warned that this may result in large computation time and memory use.
    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 122880' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 122880)' or larger];
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
    ## <td class="gt_row gt_right" style="background-color: #F1F6F0; color: #000000;">0.0065</td>
    ## <td class="gt_row gt_right">0.0036</td>
    ## <td class="gt_row gt_right">1.7785</td>
    ## <td class="gt_row gt_right">0.2836848</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #F0E7F0; color: #000000;">&minus;0.0149</td>
    ## <td class="gt_row gt_right">0.0036</td>
    ## <td class="gt_row gt_right">&minus;4.0827</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.6 &times; 10<sup class='gt_super'>&minus;4</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #BA9BCB; color: #000000;">&minus;0.0571</td>
    ## <td class="gt_row gt_right">0.0036</td>
    ## <td class="gt_row gt_right">&minus;15.6566</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">grid - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #ECDFED; color: #000000;">&minus;0.0214</td>
    ## <td class="gt_row gt_right">0.0036</td>
    ## <td class="gt_row gt_right">&minus;5.8612</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.8 &times; 10<sup class='gt_super'>&minus;8</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">grid - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #AF8DC3; color: #000000;">&minus;0.0636</td>
    ## <td class="gt_row gt_right">0.0036</td>
    ## <td class="gt_row gt_right">&minus;17.4351</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">rand - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #D5BCDC; color: #000000;">&minus;0.0422</td>
    ## <td class="gt_row gt_right">0.0036</td>
    ## <td class="gt_row gt_right">&minus;11.5739</td>
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
run_lmer(lfmm_ind, "TOTALN", filepath = here(p4path, "LFMM_individual_TOTALN.csv"))
```

    ## boundary (singular) fit: see help('isSingular')

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
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;0.1407</td>
    ## <td class="gt_row gt_right">12.2579M</td>
    ## <td class="gt_row gt_right">12.2579M</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8700K</td>
    ## <td class="gt_row gt_right">2.359458e+01</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.2 &times; 10<sup class='gt_super'>&minus;6</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">410.1322</td>
    ## <td class="gt_row gt_right">238.5254K</td>
    ## <td class="gt_row gt_right">79.5085K</td>
    ## <td class="gt_row gt_right">3</td>
    ## <td class="gt_row gt_right">122.8700K</td>
    ## <td class="gt_row gt_right">1.530411e-01</td>
    ## <td class="gt_row gt_right">0.93</td></tr>
    ##     <tr><td class="gt_row gt_left">K</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;1.2256</td>
    ## <td class="gt_row gt_right">46.1470K</td>
    ## <td class="gt_row gt_right">46.1470K</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8700K</td>
    ## <td class="gt_row gt_right">8.882559e-02</td>
    ## <td class="gt_row gt_right">0.77</td></tr>
    ##     <tr><td class="gt_row gt_left">m</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;1.9385</td>
    ## <td class="gt_row gt_right">115.4382K</td>
    ## <td class="gt_row gt_right">115.4382K</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8700K</td>
    ## <td class="gt_row gt_right">2.222002e-01</td>
    ## <td class="gt_row gt_right">0.64</td></tr>
    ##     <tr><td class="gt_row gt_left">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">1.3549</td>
    ## <td class="gt_row gt_right">56.3970K</td>
    ## <td class="gt_row gt_right">56.3970K</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8700K</td>
    ## <td class="gt_row gt_right">1.085552e-01</td>
    ## <td class="gt_row gt_right">0.74</td></tr>
    ##     <tr><td class="gt_row gt_left">H</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F7F6; color: #000000;">5.6962</td>
    ## <td class="gt_row gt_right">996.7651K</td>
    ## <td class="gt_row gt_right">996.7651K</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8700K</td>
    ## <td class="gt_row gt_right">1.918614e+00</td>
    ## <td class="gt_row gt_right">0.17</td></tr>
    ##     <tr><td class="gt_row gt_left">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">0.0219</td>
    ## <td class="gt_row gt_right">14.6781</td>
    ## <td class="gt_row gt_right">14.6781</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">122.8700K</td>
    ## <td class="gt_row gt_right">2.825307e-05</td>
    ## <td class="gt_row gt_right">1.00</td></tr>
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
    ## To enable adjustments, add the argument 'pbkrtest.limit = 122880' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 122880)' or larger];
    ## but be warned that this may result in large computation time and memory use.
    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 122880' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 122880)' or larger];
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
    ## <td class="gt_row gt_right" style="background-color: #F4F1F4; color: #000000;">&minus;0.2996</td>
    ## <td class="gt_row gt_right">5.8158</td>
    ## <td class="gt_row gt_right">&minus;0.0515</td>
    ## <td class="gt_row gt_right">0.9999509</td></tr>
    ##     <tr><td class="gt_row gt_left">envgeo - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #F2EBF2; color: #000000;">&minus;0.5890</td>
    ## <td class="gt_row gt_right">5.8158</td>
    ## <td class="gt_row gt_right">&minus;0.1013</td>
    ## <td class="gt_row gt_right">0.9996284</td></tr>
    ##     <tr><td class="gt_row gt_left">envgeo - trans</td>
    ## <td class="gt_row gt_right" style="background-color: #AF8DC3; color: #000000;">&minus;3.4776</td>
    ## <td class="gt_row gt_right">5.8158</td>
    ## <td class="gt_row gt_right">&minus;0.5980</td>
    ## <td class="gt_row gt_right">0.9327130</td></tr>
    ##     <tr><td class="gt_row gt_left">grid - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #F4F1F5; color: #000000;">&minus;0.2894</td>
    ## <td class="gt_row gt_right">5.8158</td>
    ## <td class="gt_row gt_right">&minus;0.0498</td>
    ## <td class="gt_row gt_right">0.9999558</td></tr>
    ##     <tr><td class="gt_row gt_left">grid - trans</td>
    ## <td class="gt_row gt_right" style="background-color: #B999C9; color: #000000;">&minus;3.1780</td>
    ## <td class="gt_row gt_right">5.8158</td>
    ## <td class="gt_row gt_right">&minus;0.5464</td>
    ## <td class="gt_row gt_right">0.9475315</td></tr>
    ##     <tr><td class="gt_row gt_left">rand - trans</td>
    ## <td class="gt_row gt_right" style="background-color: #C2A5D0; color: #000000;">&minus;2.8886</td>
    ## <td class="gt_row gt_right">5.8158</td>
    ## <td class="gt_row gt_right">&minus;0.4967</td>
    ## <td class="gt_row gt_right">0.9598408</td></tr>
    ##   </tbody>
    ##   
    ##   
    ## </table>
    ## </div>

### 1.1.3 Megaplots

``` r
MEGAPLOT(lfmm_ind, "K.1", colpal = "turbo")
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
MEGAPLOT(lfmm_ind, "TPRCOMBO", colpal = "plasma")
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

``` r
MEGAPLOT(lfmm_ind, "FDRCOMBO", colpal = "viridis", direction = -1)
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->

``` r
MEGAPLOT(lfmm_ind, "TOTALN", colpal = "viridis")
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-9-4.png)<!-- -->

## 1.2 Site sampling

### 1.2.1 Summary plots

### K

``` r
walk(K_plots, grid.arrange)
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-12-3.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-12-4.png)<!-- -->

### TPRCOMBO

``` r
walk(TPR_plots, grid.arrange)
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-13-4.png)<!-- -->

### FDRCOMBO

``` r
walk(FDR_plots, grid.arrange)
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-14-3.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-14-4.png)<!-- -->

### TOTALN

``` r
walk(TOTALN_plots, grid.arrange)
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->![](analysis_GEA_files/figure-gfm/unnamed-chunk-15-4.png)<!-- -->

### 1.1.2 Model summaries

``` r
lfmm_site <- 
  lfmm_site %>% 
  filter(method == "ridge") %>%
  filter(K_selection == "tracy.widom")
```

``` r
run_lmer(lfmm_site, "TPRCOMBO", filepath = here(p4path, "LFMM_site_TPR.csv"))
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
    ## <td class="gt_row gt_right" style="background-color: #F6F7F6; color: #000000;">0.0021</td>
    ## <td class="gt_row gt_right">12.4913</td>
    ## <td class="gt_row gt_right">12.4913</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">486.73649</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.7 &times; 10<sup class='gt_super'>&minus;107</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #DEF1D8; color: #000000;">0.0473</td>
    ## <td class="gt_row gt_right">2.8376</td>
    ## <td class="gt_row gt_right">1.4188</td>
    ## <td class="gt_row gt_right">2</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">55.28554</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.0 &times; 10<sup class='gt_super'>&minus;24</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #EBF4E9; color: #000000;">0.0223</td>
    ## <td class="gt_row gt_right">8.5946</td>
    ## <td class="gt_row gt_right">8.5946</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">334.89519</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.2 &times; 10<sup class='gt_super'>&minus;74</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #E2F2DE; color: #000000;">0.0386</td>
    ## <td class="gt_row gt_right">25.7362</td>
    ## <td class="gt_row gt_right">25.7362</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">1002.83784</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.6 &times; 10<sup class='gt_super'>&minus;218</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #D7EFD1; color: #000000;">0.0570</td>
    ## <td class="gt_row gt_right">56.1615</td>
    ## <td class="gt_row gt_right">56.1615</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">2188.38899</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">H</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.1112</td>
    ## <td class="gt_row gt_right">213.7780</td>
    ## <td class="gt_row gt_right">213.7780</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">8330.06915</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F0E6F0; color: #000000;">&minus;0.0262</td>
    ## <td class="gt_row gt_right">11.8362</td>
    ## <td class="gt_row gt_right">11.8362</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">461.21029</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">5.7 &times; 10<sup class='gt_super'>&minus;102</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
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
    ## To enable adjustments, add the argument 'pbkrtest.limit = 69120' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 69120)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 69120' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 69120)' or larger];
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
    ## <td class="gt_row gt_right" style="background-color: #EEE3EF; color: #000000;">&minus;0.0042</td>
    ## <td class="gt_row gt_right">0.0015</td>
    ## <td class="gt_row gt_right">&minus;2.8425</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.01244359<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - rand<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #B2DAAB; color: #000000;">0.0110</td>
    ## <td class="gt_row gt_right">0.0015</td>
    ## <td class="gt_row gt_right">7.3462</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">6.3 &times; 10<sup class='gt_super'>&minus;13</sup><sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">equi - rand<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.0152</td>
    ## <td class="gt_row gt_right">0.0015</td>
    ## <td class="gt_row gt_right">10.1887</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.7 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
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
run_lmer(lfmm_site, "FDRCOMBO", filepath = here(p4path, "LFMM_site_FDR.csv"))
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
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;0.0061</td>
    ## <td class="gt_row gt_right">109.7460</td>
    ## <td class="gt_row gt_right">109.7460</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">4782.314725</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">1.0773</td>
    ## <td class="gt_row gt_right">0.1441</td>
    ## <td class="gt_row gt_right">0.0721</td>
    ## <td class="gt_row gt_right">2</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">3.140461</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.043<sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">0.0033</td>
    ## <td class="gt_row gt_right">0.1923</td>
    ## <td class="gt_row gt_right">0.1923</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">8.377551</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.8 &times; 10<sup class='gt_super'>&minus;3</sup><sup class="gt_footnote_marks gt_asterisk">*</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">0.0043</td>
    ## <td class="gt_row gt_right">0.3241</td>
    ## <td class="gt_row gt_right">0.3241</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">14.125085</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.7 &times; 10<sup class='gt_super'>&minus;4</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F5F6; color: #000000;">&minus;0.0234</td>
    ## <td class="gt_row gt_right">9.4821</td>
    ## <td class="gt_row gt_right">9.4821</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">413.193293</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.4 &times; 10<sup class='gt_super'>&minus;91</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">H</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F6F7; color: #000000;">&minus;0.0080</td>
    ## <td class="gt_row gt_right">1.0985</td>
    ## <td class="gt_row gt_right">1.0985</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">47.867816</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">4.6 &times; 10<sup class='gt_super'>&minus;12</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">0.0034</td>
    ## <td class="gt_row gt_right">0.2004</td>
    ## <td class="gt_row gt_right">0.2004</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">8.732567</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.1 &times; 10<sup class='gt_super'>&minus;3</sup><sup class="gt_footnote_marks gt_asterisk">*</sup></td></tr>
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
    ## To enable adjustments, add the argument 'pbkrtest.limit = 69120' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 69120)' or larger];
    ## but be warned that this may result in large computation time and memory use.
    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 69120' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 69120)' or larger];
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
    ## <td class="gt_row gt_right" style="background-color: #AF8DC3; color: #000000;">&minus;0.0035</td>
    ## <td class="gt_row gt_right">0.0014</td>
    ## <td class="gt_row gt_right">&minus;2.5062</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.03269492<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">envgeo - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #E7D3E8; color: #000000;">&minus;0.0018</td>
    ## <td class="gt_row gt_right">0.0014</td>
    ## <td class="gt_row gt_right">&minus;1.2624</td>
    ## <td class="gt_row gt_right">0.41653734</td></tr>
    ##     <tr><td class="gt_row gt_left">equi - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #D9F0D3; color: #000000;">0.0018</td>
    ## <td class="gt_row gt_right">0.0014</td>
    ## <td class="gt_row gt_right">1.2437</td>
    ## <td class="gt_row gt_right">0.42738165</td></tr>
    ##   </tbody>
    ##   
    ##   <tfoot class="gt_footnotes">
    ##     <tr>
    ##       <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks gt_asterisk">***</sup> p &lt; 0.05</td>
    ##     </tr>
    ##   </tfoot>
    ## </table>
    ## </div>

``` r
run_lmer(lfmm_site, "TOTALN", filepath = here(p4path, "LFMM_site_TOTALN.csv"))
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
    ## <td class="gt_row gt_right" style="background-color: #F6F6F6; color: #000000;">&minus;16.2931</td>
    ## <td class="gt_row gt_right">786.9630M</td>
    ## <td class="gt_row gt_right">786.9630M</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">1.135564e+03</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">6.3 &times; 10<sup class='gt_super'>&minus;247</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">785.9673</td>
    ## <td class="gt_row gt_right">59.2591K</td>
    ## <td class="gt_row gt_right">29.6295K</td>
    ## <td class="gt_row gt_right">2</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">4.275454e-02</td>
    ## <td class="gt_row gt_right">0.96</td></tr>
    ##     <tr><td class="gt_row gt_left">K</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">3.3547</td>
    ## <td class="gt_row gt_right">194.4712K</td>
    ## <td class="gt_row gt_right">194.4712K</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">2.806161e-01</td>
    ## <td class="gt_row gt_right">0.60</td></tr>
    ##     <tr><td class="gt_row gt_left">m</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">0.9699</td>
    ## <td class="gt_row gt_right">16.2547K</td>
    ## <td class="gt_row gt_right">16.2547K</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">2.345501e-02</td>
    ## <td class="gt_row gt_right">0.88</td></tr>
    ##     <tr><td class="gt_row gt_left">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;1.0658</td>
    ## <td class="gt_row gt_right">19.6299K</td>
    ## <td class="gt_row gt_right">19.6299K</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">2.832532e-02</td>
    ## <td class="gt_row gt_right">0.87</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">H</td>
    ## <td class="gt_row gt_right" style="background-color: #F5F7F5; color: #000000;">20.6330</td>
    ## <td class="gt_row gt_right">7.3564M</td>
    ## <td class="gt_row gt_right">7.3564M</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">1.061510e+01</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.1 &times; 10<sup class='gt_super'>&minus;3</sup><sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;1.2650</td>
    ## <td class="gt_row gt_right">27.6501K</td>
    ## <td class="gt_row gt_right">27.6501K</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">69.1090K</td>
    ## <td class="gt_row gt_right">3.989828e-02</td>
    ## <td class="gt_row gt_right">0.84</td></tr>
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
    ## To enable adjustments, add the argument 'pbkrtest.limit = 69120' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 69120)' or larger];
    ## but be warned that this may result in large computation time and memory use.
    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 69120' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 69120)' or larger];
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
    ##     <tr><td class="gt_row gt_left">envgeo - equi</td>
    ## <td class="gt_row gt_right" style="background-color: #A4D39E; color: #000000;">1.7073</td>
    ## <td class="gt_row gt_right">7.7561</td>
    ## <td class="gt_row gt_right">0.2201</td>
    ## <td class="gt_row gt_right">0.9736425</td></tr>
    ##     <tr><td class="gt_row gt_left">envgeo - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">2.1467</td>
    ## <td class="gt_row gt_right">7.7561</td>
    ## <td class="gt_row gt_right">0.2768</td>
    ## <td class="gt_row gt_right">0.9586535</td></tr>
    ##     <tr><td class="gt_row gt_left">equi - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #EBF4E8; color: #000000;">0.4394</td>
    ## <td class="gt_row gt_right">7.7561</td>
    ## <td class="gt_row gt_right">0.0566</td>
    ## <td class="gt_row gt_right">0.9982324</td></tr>
    ##   </tbody>
    ##   
    ##   
    ## </table>
    ## </div>

### 1.2.3 Megaplots

``` r
MEGAPLOT(lfmm_site, "K.1", colpal = "turbo")
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
MEGAPLOT(lfmm_site, "TPRCOMBO", colpal = "plasma")
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

``` r
MEGAPLOT(lfmm_site, "FDRCOMBO", colpal = "viridis", direction = -1)
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-18-3.png)<!-- -->

``` r
MEGAPLOT(lfmm_site, "TOTALN", colpal = "viridis")
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-18-4.png)<!-- -->

# 2. RDA

## 2.1 Individual sampling

### 2.1.1 Summary plots

``` r
rda_ind <- format_rda(here(dirname(getwd()), "p3_methods", "outputs", "rda_indsampling_results.csv"))

summary_hplot(rda_ind, "TPRCOMBO", colpal = "plasma")
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
summary_hplot(rda_ind, "FDRCOMBO", colpal = "viridis", direction = -1)
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

``` r
summary_hplot(rda_ind, "TOTALN", colpal = "viridis", direction = -1)
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-19-3.png)<!-- -->

### 2.1.2 Model summaries

``` r
run_lmer(rda_ind, "TPRCOMBO", filepath = here(p4path, "RDA_individual_TPR.csv"))
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
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">0.0002</td>
    ## <td class="gt_row gt_right">5.9014</td>
    ## <td class="gt_row gt_right">5.9014</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">1046.25067</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">9.8 &times; 10<sup class='gt_super'>&minus;226</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.0915</td>
    ## <td class="gt_row gt_right">0.6548</td>
    ## <td class="gt_row gt_right">0.2183</td>
    ## <td class="gt_row gt_right">3</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">38.69771</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">6.0 &times; 10<sup class='gt_super'>&minus;25</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #ECF5EA; color: #000000;">0.0167</td>
    ## <td class="gt_row gt_right">2.1375</td>
    ## <td class="gt_row gt_right">2.1375</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">378.95346</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">6.7 &times; 10<sup class='gt_super'>&minus;84</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #D7EFD1; color: #000000;">0.0469</td>
    ## <td class="gt_row gt_right">16.8867</td>
    ## <td class="gt_row gt_right">16.8867</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">2993.81298</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #D6EFD0; color: #000000;">0.0471</td>
    ## <td class="gt_row gt_right">17.0395</td>
    ## <td class="gt_row gt_right">17.0395</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">3020.89204</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">H</td>
    ## <td class="gt_row gt_right" style="background-color: #D7EFD1; color: #000000;">0.0467</td>
    ## <td class="gt_row gt_right">16.7580</td>
    ## <td class="gt_row gt_right">16.7580</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">2970.99516</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F5F6F4; color: #000000;">0.0035</td>
    ## <td class="gt_row gt_right">0.0949</td>
    ## <td class="gt_row gt_right">0.0949</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">16.82851</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">4.1 &times; 10<sup class='gt_super'>&minus;5</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
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
    ## To enable adjustments, add the argument 'pbkrtest.limit = 30720' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 30720)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 30720' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 30720)' or larger];
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
    ## <td class="gt_row gt_right" style="background-color: #F5F7F5; color: #000000;">0.0003</td>
    ## <td class="gt_row gt_right">0.0012</td>
    ## <td class="gt_row gt_right">0.2552</td>
    ## <td class="gt_row gt_right">0.99417692</td></tr>
    ##     <tr><td class="gt_row gt_left">envgeo - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #E7F3E4; color: #000000;">0.0031</td>
    ## <td class="gt_row gt_right">0.0012</td>
    ## <td class="gt_row gt_right">2.5381</td>
    ## <td class="gt_row gt_right">0.05424359</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.0114</td>
    ## <td class="gt_row gt_right">0.0012</td>
    ## <td class="gt_row gt_right">9.4274</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.6 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">grid - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #E9F4E6; color: #000000;">0.0028</td>
    ## <td class="gt_row gt_right">0.0012</td>
    ## <td class="gt_row gt_right">2.2830</td>
    ## <td class="gt_row gt_right">0.10196730</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">grid - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #84C280; color: #000000;">0.0111</td>
    ## <td class="gt_row gt_right">0.0012</td>
    ## <td class="gt_row gt_right">9.1722</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">4.7 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">rand - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #B0D9AA; color: #000000;">0.0083</td>
    ## <td class="gt_row gt_right">0.0012</td>
    ## <td class="gt_row gt_right">6.8892</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.4 &times; 10<sup class='gt_super'>&minus;11</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
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
run_lmer(rda_ind, "FDRCOMBO", filepath = here(p4path, "RDA_individual_FDR.csv"))
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
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">0.0003</td>
    ## <td class="gt_row gt_right">12.4744</td>
    ## <td class="gt_row gt_right">12.4744</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">671.624779</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.7 &times; 10<sup class='gt_super'>&minus;146</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.1207</td>
    ## <td class="gt_row gt_right">0.5150</td>
    ## <td class="gt_row gt_right">0.1717</td>
    ## <td class="gt_row gt_right">3</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">9.241763</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">4.2 &times; 10<sup class='gt_super'>&minus;6</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #F0F5EE; color: #000000;">0.0150</td>
    ## <td class="gt_row gt_right">1.7308</td>
    ## <td class="gt_row gt_right">1.7308</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">93.187692</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">5.1 &times; 10<sup class='gt_super'>&minus;22</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #B4DCAE; color: #000000;">0.0852</td>
    ## <td class="gt_row gt_right">55.7350</td>
    ## <td class="gt_row gt_right">55.7350</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">3000.793965</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #CDEAC7; color: #000000;">0.0683</td>
    ## <td class="gt_row gt_right">35.7750</td>
    ## <td class="gt_row gt_right">35.7750</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">1926.140142</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">H</td>
    ## <td class="gt_row gt_right" style="background-color: #D3EDCD; color: #000000;">0.0644</td>
    ## <td class="gt_row gt_right">31.8078</td>
    ## <td class="gt_row gt_right">31.8078</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">1712.547206</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F3F6F3; color: #000000;">0.0073</td>
    ## <td class="gt_row gt_right">0.4050</td>
    ## <td class="gt_row gt_right">0.4050</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">21.804129</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.0 &times; 10<sup class='gt_super'>&minus;6</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
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
    ## To enable adjustments, add the argument 'pbkrtest.limit = 30720' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 30720)' or larger];
    ## but be warned that this may result in large computation time and memory use.
    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 30720' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 30720)' or larger];
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
    ## <td class="gt_row gt_right" style="background-color: #F3EFF3; color: #000000;">&minus;0.0012</td>
    ## <td class="gt_row gt_right">0.0022</td>
    ## <td class="gt_row gt_right">&minus;0.5647</td>
    ## <td class="gt_row gt_right">0.942526605</td></tr>
    ##     <tr><td class="gt_row gt_left">envgeo - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #E9F4E6; color: #000000;">0.0025</td>
    ## <td class="gt_row gt_right">0.0022</td>
    ## <td class="gt_row gt_right">1.1249</td>
    ## <td class="gt_row gt_right">0.674115705</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #94CB8F; color: #000000;">0.0093</td>
    ## <td class="gt_row gt_right">0.0022</td>
    ## <td class="gt_row gt_right">4.2501</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.3 &times; 10<sup class='gt_super'>&minus;4</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">grid - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #E2F2DE; color: #000000;">0.0037</td>
    ## <td class="gt_row gt_right">0.0022</td>
    ## <td class="gt_row gt_right">1.6896</td>
    ## <td class="gt_row gt_right">0.329123417</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">grid - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.0106</td>
    ## <td class="gt_row gt_right">0.0022</td>
    ## <td class="gt_row gt_right">4.8148</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">8.8 &times; 10<sup class='gt_super'>&minus;6</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">rand - trans<sup class="gt_footnote_marks gt_asterisk">**</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #BEE1B8; color: #000000;">0.0069</td>
    ## <td class="gt_row gt_right">0.0022</td>
    ## <td class="gt_row gt_right">3.1252</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.009619135<sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
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
run_lmer(rda_ind, "TOTALN", filepath = here(p4path, "RDA_individual_TOTALN.csv"))
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
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">0.0030</td>
    ## <td class="gt_row gt_right">1.4041K</td>
    ## <td class="gt_row gt_right">1.4041K</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">1070.36871</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">8.4 &times; 10<sup class='gt_super'>&minus;231</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">1.3736</td>
    ## <td class="gt_row gt_right">131.2011</td>
    ## <td class="gt_row gt_right">43.7337</td>
    ## <td class="gt_row gt_right">3</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">33.33937</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.7 &times; 10<sup class='gt_super'>&minus;21</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #ECF5EA; color: #000000;">0.2472</td>
    ## <td class="gt_row gt_right">469.3102</td>
    ## <td class="gt_row gt_right">469.3102</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">357.76763</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.4 &times; 10<sup class='gt_super'>&minus;79</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #D5EECF; color: #000000;">0.7201</td>
    ## <td class="gt_row gt_right">3.9826K</td>
    ## <td class="gt_row gt_right">3.9826K</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">3036.04799</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #D9F0D3; color: #000000;">0.6906</td>
    ## <td class="gt_row gt_right">3.6624K</td>
    ## <td class="gt_row gt_right">3.6624K</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">2791.93296</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">H</td>
    ## <td class="gt_row gt_right" style="background-color: #D9F0D3; color: #000000;">0.6827</td>
    ## <td class="gt_row gt_right">3.5800K</td>
    ## <td class="gt_row gt_right">3.5800K</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">2729.11844</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F3F6F2; color: #000000;">0.0890</td>
    ## <td class="gt_row gt_right">60.8297</td>
    ## <td class="gt_row gt_right">60.8297</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">30.7080K</td>
    ## <td class="gt_row gt_right">46.37211</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.0 &times; 10<sup class='gt_super'>&minus;11</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
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
    ## To enable adjustments, add the argument 'pbkrtest.limit = 30720' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 30720)' or larger];
    ## but be warned that this may result in large computation time and memory use.
    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 30720' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 30720)' or larger];
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
    ## <td class="gt_row gt_right" style="background-color: #F5F7F5; color: #000000;">0.0046</td>
    ## <td class="gt_row gt_right">0.0185</td>
    ## <td class="gt_row gt_right">0.2466</td>
    ## <td class="gt_row gt_right">0.99473672</td></tr>
    ##     <tr><td class="gt_row gt_left">envgeo - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #E6F3E3; color: #000000;">0.0458</td>
    ## <td class="gt_row gt_right">0.0185</td>
    ## <td class="gt_row gt_right">2.4798</td>
    ## <td class="gt_row gt_right">0.06307795</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.1620</td>
    ## <td class="gt_row gt_right">0.0185</td>
    ## <td class="gt_row gt_right">8.7639</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.1 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left">grid - rand</td>
    ## <td class="gt_row gt_right" style="background-color: #E8F4E5; color: #000000;">0.0413</td>
    ## <td class="gt_row gt_right">0.0185</td>
    ## <td class="gt_row gt_right">2.2332</td>
    ## <td class="gt_row gt_right">0.11431421</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">grid - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #84C280; color: #000000;">0.1574</td>
    ## <td class="gt_row gt_right">0.0185</td>
    ## <td class="gt_row gt_right">8.5173</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.7 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">rand - trans<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #B2DBAC; color: #000000;">0.1161</td>
    ## <td class="gt_row gt_right">0.0185</td>
    ## <td class="gt_row gt_right">6.2841</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.0 &times; 10<sup class='gt_super'>&minus;9</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
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
MEGAPLOT(rda_ind, "TPRCOMBO", colpal = "plasma")
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
MEGAPLOT(rda_ind, "FDRCOMBO", colpal = "viridis", direction = -1)
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-21-2.png)<!-- -->

``` r
MEGAPLOT(rda_ind, "TOTALN", colpal = "viridis")
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-21-3.png)<!-- -->

## 2.2 Site sampling

### 2.2.1 Summary plots

``` r
rda_site <- format_rda(here(dirname(getwd()), "p3_methods", "outputs", "rda_sitesampling_results.csv"))

summary_hplot(rda_site, "TPRCOMBO", colpal = "plasma")
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
summary_hplot(rda_site, "FDRCOMBO", colpal = "viridis", direction = -1)
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-22-2.png)<!-- -->

``` r
summary_hplot(rda_site, "TOTALN", colpal = "viridis")
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-22-3.png)<!-- -->

### 2.2.2 Model summaries

``` r
run_lmer(rda_site, "TPRCOMBO", filepath = here(p4path, "RDA_site_TPR.csv"))
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
    ## <td class="gt_row gt_right" style="background-color: #F6F7F6; color: #000000;">0.0016</td>
    ## <td class="gt_row gt_right">1.8569</td>
    ## <td class="gt_row gt_right">1.8569</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">431.692944</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.0 &times; 10<sup class='gt_super'>&minus;94</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.0798</td>
    ## <td class="gt_row gt_right">0.2606</td>
    ## <td class="gt_row gt_right">0.1303</td>
    ## <td class="gt_row gt_right">2</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">30.287115</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">7.4 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #E9F4E5; color: #000000;">0.0193</td>
    ## <td class="gt_row gt_right">1.6164</td>
    ## <td class="gt_row gt_right">1.6164</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">375.762993</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">7.9 &times; 10<sup class='gt_super'>&minus;83</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #DBF1D5; color: #000000;">0.0372</td>
    ## <td class="gt_row gt_right">5.9677</td>
    ## <td class="gt_row gt_right">5.9677</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">1387.328794</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.8 &times; 10<sup class='gt_super'>&minus;292</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #DBF1D6; color: #000000;">0.0371</td>
    ## <td class="gt_row gt_right">5.9306</td>
    ## <td class="gt_row gt_right">5.9306</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">1378.701819</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.1 &times; 10<sup class='gt_super'>&minus;290</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">H</td>
    ## <td class="gt_row gt_right" style="background-color: #DCF1D6; color: #000000;">0.0364</td>
    ## <td class="gt_row gt_right">5.7376</td>
    ## <td class="gt_row gt_right">5.7376</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">1333.851639</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.3 &times; 10<sup class='gt_super'>&minus;281</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F5F7F5; color: #000000;">0.0027</td>
    ## <td class="gt_row gt_right">0.0323</td>
    ## <td class="gt_row gt_right">0.0323</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">7.508874</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">6.1 &times; 10<sup class='gt_super'>&minus;3</sup><sup class="gt_footnote_marks gt_asterisk">**</sup></td></tr>
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
    ## To enable adjustments, add the argument 'pbkrtest.limit = 17280' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 17280)' or larger];
    ## but be warned that this may result in large computation time and memory use.

    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 17280' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 17280)' or larger];
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
    ##     <tr><td class="gt_row gt_left">envgeo - equi</td>
    ## <td class="gt_row gt_right" style="background-color: #F6F5F6; color: #000000;">&minus;0.0002</td>
    ## <td class="gt_row gt_right">0.0012</td>
    ## <td class="gt_row gt_right">&minus;0.1598</td>
    ## <td class="gt_row gt_right">0.9860183</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #83C17F; color: #000000;">0.0081</td>
    ## <td class="gt_row gt_right">0.0012</td>
    ## <td class="gt_row gt_right">6.6589</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">8.3 &times; 10<sup class='gt_super'>&minus;11</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">equi - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.0083</td>
    ## <td class="gt_row gt_right">0.0012</td>
    ## <td class="gt_row gt_right">6.8187</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.8 &times; 10<sup class='gt_super'>&minus;11</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
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
run_lmer(rda_site, "FDRCOMBO", filepath = here(p4path, "RDA_site_FDR.csv"))
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
    ## <td class="gt_row gt_right" style="background-color: #F6F7F6; color: #000000;">0.0025</td>
    ## <td class="gt_row gt_right">4.7325</td>
    ## <td class="gt_row gt_right">4.7325</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">330.96656</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.9 &times; 10<sup class='gt_super'>&minus;73</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.1289</td>
    ## <td class="gt_row gt_right">0.5401</td>
    ## <td class="gt_row gt_right">0.2700</td>
    ## <td class="gt_row gt_right">2</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">18.88484</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">6.4 &times; 10<sup class='gt_super'>&minus;9</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #EDF5EA; color: #000000;">0.0225</td>
    ## <td class="gt_row gt_right">2.1796</td>
    ## <td class="gt_row gt_right">2.1796</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">152.43043</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">7.2 &times; 10<sup class='gt_super'>&minus;35</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #CEEAC8; color: #000000;">0.0726</td>
    ## <td class="gt_row gt_right">22.7847</td>
    ## <td class="gt_row gt_right">22.7847</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">1593.45311</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">0.0<sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #D7EFD1; color: #000000;">0.0659</td>
    ## <td class="gt_row gt_right">18.7672</td>
    ## <td class="gt_row gt_right">18.7672</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">1312.48735</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">4.6 &times; 10<sup class='gt_super'>&minus;277</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">H</td>
    ## <td class="gt_row gt_right" style="background-color: #DBF0D5; color: #000000;">0.0602</td>
    ## <td class="gt_row gt_right">15.6564</td>
    ## <td class="gt_row gt_right">15.6564</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">1094.93525</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">7.4 &times; 10<sup class='gt_super'>&minus;233</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F3F6F2; color: #000000;">0.0095</td>
    ## <td class="gt_row gt_right">0.3888</td>
    ## <td class="gt_row gt_right">0.3888</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">27.18806</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.9 &times; 10<sup class='gt_super'>&minus;7</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
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
    ## To enable adjustments, add the argument 'pbkrtest.limit = 17280' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 17280)' or larger];
    ## but be warned that this may result in large computation time and memory use.
    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 17280' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 17280)' or larger];
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
    ##     <tr><td class="gt_row gt_left">envgeo - equi</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F7F7; color: #000000;">&minus;0.0001</td>
    ## <td class="gt_row gt_right">0.0022</td>
    ## <td class="gt_row gt_right">&minus;0.0319</td>
    ## <td class="gt_row gt_right">0.9994395</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #80C07C; color: #000000;">0.0118</td>
    ## <td class="gt_row gt_right">0.0022</td>
    ## <td class="gt_row gt_right">5.3063</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.4 &times; 10<sup class='gt_super'>&minus;7</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">equi - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.0119</td>
    ## <td class="gt_row gt_right">0.0022</td>
    ## <td class="gt_row gt_right">5.3382</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.8 &times; 10<sup class='gt_super'>&minus;7</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
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
run_lmer(rda_site, "TOTALN", filepath = here(p4path, "RDA_site_TOTALN.csv"))
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
    ## <td class="gt_row gt_right" style="background-color: #F6F7F6; color: #000000;">0.0247</td>
    ## <td class="gt_row gt_right">451.8537</td>
    ## <td class="gt_row gt_right">451.8537</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">438.65917</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.3 &times; 10<sup class='gt_super'>&minus;96</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">sampstrat</td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">1.2344</td>
    ## <td class="gt_row gt_right">65.8091</td>
    ## <td class="gt_row gt_right">32.9046</td>
    ## <td class="gt_row gt_right">2</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">31.94373</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.4 &times; 10<sup class='gt_super'>&minus;14</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">K</td>
    ## <td class="gt_row gt_right" style="background-color: #E9F4E6; color: #000000;">0.2932</td>
    ## <td class="gt_row gt_right">371.3014</td>
    ## <td class="gt_row gt_right">371.3014</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">360.45913</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.4 &times; 10<sup class='gt_super'>&minus;79</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">m</td>
    ## <td class="gt_row gt_right" style="background-color: #DBF0D5; color: #000000;">0.5760</td>
    ## <td class="gt_row gt_right">1.4335K</td>
    ## <td class="gt_row gt_right">1.4335K</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">1391.62088</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">5.2 &times; 10<sup class='gt_super'>&minus;293</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">phi</td>
    ## <td class="gt_row gt_right" style="background-color: #DCF1D6; color: #000000;">0.5626</td>
    ## <td class="gt_row gt_right">1.3674K</td>
    ## <td class="gt_row gt_right">1.3674K</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">1327.50723</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">4.3 &times; 10<sup class='gt_super'>&minus;280</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">H</td>
    ## <td class="gt_row gt_right" style="background-color: #DCF1D7; color: #000000;">0.5478</td>
    ## <td class="gt_row gt_right">1.2964K</td>
    ## <td class="gt_row gt_right">1.2964K</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">1258.51577</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">3.8 &times; 10<sup class='gt_super'>&minus;266</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">r</td>
    ## <td class="gt_row gt_right" style="background-color: #F4F6F3; color: #000000;">0.0652</td>
    ## <td class="gt_row gt_right">18.3431</td>
    ## <td class="gt_row gt_right">18.3431</td>
    ## <td class="gt_row gt_right">1</td>
    ## <td class="gt_row gt_right">17.2690K</td>
    ## <td class="gt_row gt_right">17.80748</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">2.5 &times; 10<sup class='gt_super'>&minus;5</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
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
    ## To enable adjustments, add the argument 'pbkrtest.limit = 17280' (or larger)
    ## [or, globally, 'set emm_options(pbkrtest.limit = 17280)' or larger];
    ## but be warned that this may result in large computation time and memory use.
    ## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
    ## To enable adjustments, add the argument 'lmerTest.limit = 17280' (or larger)
    ## [or, globally, 'set emm_options(lmerTest.limit = 17280)' or larger];
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
    ##     <tr><td class="gt_row gt_left">envgeo - equi</td>
    ## <td class="gt_row gt_right" style="background-color: #F7F6F7; color: #000000;">&minus;0.0017</td>
    ## <td class="gt_row gt_right">0.0189</td>
    ## <td class="gt_row gt_right">&minus;0.0918</td>
    ## <td class="gt_row gt_right">0.9953648</td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">envgeo - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #81C07D; color: #000000;">0.1300</td>
    ## <td class="gt_row gt_right">0.0189</td>
    ## <td class="gt_row gt_right">6.8758</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">1.9 &times; 10<sup class='gt_super'>&minus;11</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
    ##     <tr><td class="gt_row gt_left" style="font-weight: bold;">equi - rand<sup class="gt_footnote_marks gt_asterisk">***</sup></td>
    ## <td class="gt_row gt_right" style="background-color: #7FBF7B; color: #000000;">0.1318</td>
    ## <td class="gt_row gt_right">0.0189</td>
    ## <td class="gt_row gt_right">6.9676</td>
    ## <td class="gt_row gt_right" style="font-weight: bold;">9.7 &times; 10<sup class='gt_super'>&minus;12</sup><sup class="gt_footnote_marks gt_asterisk">***</sup></td></tr>
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
MEGAPLOT(rda_site, "TPRCOMBO", colpal = "plasma")
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
MEGAPLOT(rda_site, "FDRCOMBO", colpal = "viridis", direction = -1)
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-24-2.png)<!-- -->

``` r
MEGAPLOT(rda_site, "TOTALN", colpal = "viridis")
```

![](analysis_GEA_files/figure-gfm/unnamed-chunk-24-3.png)<!-- -->
