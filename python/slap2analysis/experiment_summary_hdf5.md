```
📦experiment_summary.h5
 ┗ 📂/
   ┣ 📂DMD1/
   ┃ ┣ 📂frame_info/
   ┃ ┃ ┣ 📜offlineXshifts (total frames x 1)
   ┃ ┃ ┣ 📜offlineYshifts (total frames x 1)
   ┃ ┃ ┣ 📜offlineZshifts (total frames x 1)
   ┃ ┃ ┣ 📜onlineXshifts (total frames x 1)
   ┃ ┃ ┣ 📜onlineYshifts (total frames x 1)
   ┃ ┃ ┣ 📜onlineZshifts (total frames x 1)
   ┃ ┃ ┣ 📜trial_num_frames (trials x 1)
   ┃ ┃ ┗ 📜discard_frames (total frames x 1)
   ┃ ┣ 📂visualizations/
   ┃ ┃ ┣ 📜mean_im (channels x fastz x rows x cols)
   ┃ ┃ ┣ 📜ref_stack (ref_stack_channels x depths x rows x cols)
   ┃ ┃ ┃  ┗ ℹ️channels (ref_stack_channels x 1)
   ┃ ┃ ┗ 📜act_im (fastz x rows x cols)
   ┃ ┣ 📂global/
   ┃ ┃ ┗ 📜F(total frames x 1)
   ┃ ┣ 📂user_rois/
   ┃ ┃ ┣ 📜mask (fastz x rows x cols x rois)
   ┃ ┃ ┣ 📜Fsvd (total frames x rois)
   ┃ ┃ ┗ 📜F (total frames x rois)
   ┃ ┗ 📂sources/
   ┃   ┣ 📂spatial/
   ┃   ┃ ┣ 📜profiles (sources x fastz x rows x cols)
   ┃   ┃ ┗ 📜coords (sources x 2 [z_loc, x_loc, y_loc])*
   ┃   ┗ 📂temporal/
   ┃     ┣ 📜dF (total frames x sources)
   ┃     ┣ 📜dFF (total frames x sources)
   ┃     ┗ 📜F0 (total frames x sources)
   ┣ 📂DMD2/
   ┃ ┣ ... (same as DMD1)
```
