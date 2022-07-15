# lv_seg
Automatic segmentation of the left ventricle endocardium on each diastolic and systolic image from short axis cardiac magnetic resonance images

Diagnosing heart diseases from imaging techniques is one of many challenges faced by
cardiologists.Segmentation of the left ventricle from Cardiac Magnetic Resonance Imaging (CMR) images
is a critical step in this process. In this report, we study and implement an automatic and unsupervised
segmentation technique for the left ventricle endocardium to achieve fast and reliable performance in
clinical biomedical applications based on the work of Yang et al. This article presents a fully automated
left ventricular segmentation in short-axis magnetic resonance imaging. Segmentation was made possible
by the use of two clustering-based algorithms; a standard Fuzzy C-Means (FCM) and Circular Shape-
constrained FCM (CS-FCM). The algorithm along with the mentioned clustering techniques has been shown
to be capable of left ventricle segmentation in the systolic and diastolic phases of the heartbeat cycle. We
have been able to reproduce said algorithm on image data obtained from “Automatic Cardiac Diagnosis
Challenge” (ACDC) data-set. Furthermore, we have also been able to improve the LV detection by means of
incorporating a hyper-parameter tuning operation.Comparisons of our results with ground truth data yielded
a 0.89 Dice similarity suggesting the method is able to effectively segment, in most cases, the LV-region
both in diastole and systole cycles.

![image](https://user-images.githubusercontent.com/19288227/179282607-8637a1d1-2ce3-415a-8edb-7aa23ea4a5fd.png)

![image](https://user-images.githubusercontent.com/19288227/179282632-102f1a84-1127-43d6-ad0d-165fd5fc4e2e.png)

![image](https://user-images.githubusercontent.com/19288227/179282678-23066609-9395-43f2-a5fa-07967a131b23.png)
