```
📦experiment_summary.h5
 ┗ 📂/
   ┣ 📂DMD1/
   ┃ ┣ 📂frame_info/
   ┃ ┃ ┣ 📜*offlineXshift (total frames x 1)*
   ┃ ┃ ┣ 📜*offlineYshift (total frames x 1)*
   ┃ ┃ ┣ 📜*offlineZshift (total frames x 1)*
   ┃ ┃ ┣ 📜*onlineXshift (total frames x 1)*
   ┃ ┃ ┣ 📜*onlineYshift (total frames x 1)*
   ┃ ┃ ┣ 📜*onlineZshift (total frames x 1)*
   ┃ ┃ ┣ 📜trial_start_idxs (trials x 1)
   ┃ ┃ ┗ 📜discard_frames (total frames x 1)
   ┃ ┣ 📂visualizations/
   ┃ ┃ ┣ 📜mean_im (rows x cols)
   ┃ ┃ ┣ 📜act_im (rows x cols)
   ┃ ┃ ┣ 📜per_trial_mean_im (trials x rows x cols)
   ┃ ┃ ┗ 📜per_trial_act_im (trials x rows x cols)
   ┃ ┣ 📂global/
   ┃ ┃ ┗ 📜F(total frames x 1)
   ┃ ┣ 📂user_rois/
   ┃ ┃ ┣ 📜mask (rows x cols x rois)
   ┃ ┃ ┣ 📜Fsvd (total frames x rois)
   ┃ ┃ ┗ 📜F (total frames x rois)
   ┃ ┗ 📂sources/
   ┃   ┣ 📂spatial/
   ┃   ┃ ┣ 📜footprints (pixels x sources)
   ┃   ┃ ┗ 📜*source_params (sources x 4 [x_loc, y_loc, x_sigma, y_sigma])*
   ┃   ┗ 📂temporal/
   ┃     ┣ 📜dF (total frames x sources)
   ┃     ┣ 📜dFF (total frames x sources)
   ┃     ┗ 📜F0 (total frames x sources)
   ┣ 📂DMD2/
   ┃ ┣ ... (same as DMD1)
```
