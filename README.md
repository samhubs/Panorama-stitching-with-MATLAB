
<p>The submitted routines can be used to stitch panoramas using any number of input images. There are two main files depending on the input images if colored or gray-scale-- "SIFT_allImages_stitch_gray" and "SIFT_allImages_stitch_new_colored".
Description along with the procedure :
1) Detect features : The most stable Harris points estimated using the Scale-space maximization of the Laplacian of Gaussian (LOG) formulation. The estimation is carried out using the function "Harris_Laplace_fn.m"
2) Using VLFEAT associate SIFT descriptor to the identified corners.
3) Match the corresponding points using the SIFT descriptors.
4) Compute Homography using the matched points. The function is "RANSAC.m"
5) Stitch images using the CANVAS approach
- > The figures are defined at specific locations and can be enabled to visualize the corners or the matched points.<p>


<div class="content"><h2>Contents</h2><div><ul><li><a href="#2">This function evaluates the most stable Harris corner points at different scales based on Laplacian of Gaussian (LOG) formulation using Scale pyramid</a></li><li><a href="#3">References:</a></li><li><a href="#4"><b>Inputs:</b></a></li><li><a href="#5"><b>Outputs:</b></a></li><li><a href="#7">This is the main loop which runs over all the images.</a></li><li><a href="#9">Part #1 , Finding Harris-corner points by thresholding <b>Harris measure,R</b></a></li><li><a href="#30">Part #2 , Checking the stability of intereset points using Laplacian of Gaussian (LOG)</a></li></ul></div><pre class="codeinput"><span class="keyword">function</span> my_final = Harris_Laplace_fn(Img, threshold_R)
</pre><h2>This function evaluates the most stable Harris corner points at different scales based on Laplacian of Gaussian (LOG) formulation using Scale pyramid<a name="2"></a></h2><p>Written by : Sameer Sharma (ED13D007), Engg. design , IITM</p><h2>References:<a name="3"></a></h2><div><ol><li>Harris-Affine &amp; Hessian Affine: K. Mikolajczyk and C. Schmid, Scale and Affine invariant interest point detectors. In IJC V 60(1):63-86, 2004.</li><li>Digital Video Processing (<a href="http://10.6.4.152/dvp/dvp.html">http://10.6.4.152/dvp/dvp.html</a>), Computer Sc., Indian Institute of Technology, Madras.</li><li>A Comparison of Affine Region Detectors, K. Mikolajczyk, T. Tuytelaars, C. Schmid, A. Zisserman, J. Matas, F. Schaffalitzky, T. Kadir, L. Van Gool. IJCV 2005.</li><li>Scale-Space Theory in Computer Vision, T. Lindeberg. Springer</li><li><a href="http://en.wikipedia.org/wiki/Corner_detection">http://en.wikipedia.org/wiki/Corner_detection</a> The_Harris_Stephens_Plessey_Shi_Tomasi_corner_detection_algorithm</li><li>Local Invariant Feature Detectors - A Survey, T. Tuytelaars and K. Mikolajczyk. FGV 2008</li></ol></div><h2><b>Inputs:</b><a name="4"></a></h2><pre>        *Img* : 'String datatype of the current image name ( If colored,
        convert to gray scale first)
        *threshold_R* : Threshold on the Harris measure for corner
        detection ( Det(M) - k* (trace(M))^2 )</pre><h2><b>Outputs:</b><a name="5"></a></h2><pre>         *my_final* : a N x 3 array with spatial cordinates and scales of
         the computed Harris interest points</pre><p>Defining scale-levels and initial-scale(these can be further customized)</p><pre class="codeinput">levels = 12;
sigma_initial = 1.2;
scale_factor = 1.2;
myPoints = [];
my_final = [];
</pre><h2>This is the main loop which runs over all the images.<a name="7"></a></h2><p>The loop is divided into segments such that it will be easy for the user to track the process.</p><pre class="codeinput"><span class="keyword">for</span> n = 1:levels
</pre><p>Defining the current scale ( Sigma value for the Guassian mask)</p><pre class="codeinput">sigma = (scale_factor)^n*sigma_initial;
</pre><h2>Part #1 , Finding Harris-corner points by thresholding <b>Harris measure,R</b><a name="9"></a></h2><p>Define a Guassian filter mask for a given value of sigma scales Making sure that the filter size is ODD with pixel of interest at the center</p><pre class="codeinput">patch_size = 6*ceil(sigma)+1;
<span class="keyword">if</span> mod(patch_size,2) ==0;
    patch_size = patch_size +1;
<span class="keyword">end</span>
</pre><p>Mask definition</p><pre class="codeinput">mask = fspecial(<span class="string">'gaussian'</span>,patch_size,sigma);
</pre><p>Gaussian differentitation (0.7 * sigma)</p><pre class="codeinput">G1 = fspecial(<span class="string">'gaussian'</span>,patch_size,0.7*sigma);
[Gx,Gy] = gradient(G1);
</pre><p>Computing second order derivatives for later use( Laplacian computation)</p><pre class="codeinput">[Gxx,~] = gradient(Gx);
[~,Gyy] = gradient(Gy);
myImage1 = double((imread(Img)));
myImage = filter2(mask, myImage1);
</pre><p>Compute the gradients of the Gaussians</p><pre class="codeinput">[ro, col] = size(myImage);
grad_x = filter2(Gx, myImage1);
grad_y = filter2(Gy, myImage1);
grad_xx = filter2(Gxx, myImage1);
grad_yy = filter2(Gyy, myImage1);
</pre><p>Assembling the Moment matrix (2 X 2)</p><pre class="codeinput">M1 = grad_x.^2;
M2 = grad_y.^2;
M3 = grad_x.*grad_y;
</pre><p>Mask the three matrices using Guassian window</p><pre class="codeinput">M1_new = filter2(mask, M1);
M2_new = filter2(mask, M2);
M3_new = filter2(mask, M3);
R = zeros(ro, col);
</pre><p>Computing the equivalent to the eigen values of the moment matrix for each pixel by defining the measure R <i><b>(Szieliski) : det(M)/trace(M)</b></i> or <i>*Harris = det(M) - K*trace(M)^2</i>*, K can be taken as 0.7</p><pre class="codeinput">R = (M1_new.*M2_new - M3_new.^2)./(M1_new + M2_new);
</pre><p>Normalize the R matrix such that the values lies [ 0 % , 100 %]</p><pre class="codeinput">R_norm = 100*(R - min(min(R)))/(max(max(R)) - min(min(R)));
</pre><p>Thresholding based on the normalised R value</p><pre class="codeinput">myLogic = R_norm &gt; threshold_R;
<span class="comment">% figure, imshow(myLogic);</span>
</pre><p>Padding the boudaries</p><pre class="codeinput">myLogic(1: ceil(patch_size/2),:) = 0;
myLogic(:,1: ceil(patch_size/2)) = 0;
myLogic(end-ceil(patch_size/2):end,:) = 0;
myLogic(:,end- ceil(patch_size/2):end) = 0;
R_norm_new = R_norm.*double(myLogic);
bound = floor(patch_size/2);
</pre><p>Check the <i>interest point</i>, if it has attained maximum value of R as compared to the neigbours i_z and index are padded to avoid boundary problems</p><pre class="codeinput"><span class="keyword">for</span> i_z = bound + 1:ro - bound
</pre><p>Scan row-wise of the thresholded masked Non-Zero R values</p><pre class="codeinput">    index = find(R_norm_new(i_z,:));
</pre><p>Check if we have an interset point at that index</p><pre class="codeinput">    size_new = size(index,2);
</pre><p>Check if the intereset point has maximum R value as compared to the neighbours</p><pre class="codeinput">    <span class="keyword">if</span> ~isempty(index) &amp;&amp; min(index)&gt; bound &amp;&amp; max(index)&lt;col - bound

        <span class="keyword">for</span> check_index = 1:size_new
</pre><pre class="codeinput">        neigbors = R_norm_new(i_z - bound:i_z +bound,<span class="keyword">...</span>
               index(check_index)- bound:index(check_index)+<span class="keyword">...</span>
              bound);
         [~, colIdx] = max(max(neigbors,[],1));
         [~, rowIdx] = max(max(neigbors,[],2));
</pre><p>Update the final mask if it is maximum</p><pre class="codeinput">            <span class="keyword">if</span> rowIdx*colIdx ~= (bound + 1)^2

            myLogic(i_z,index(check_index)) = 0;
            R_norm_new(i_z,index(check_index)) = 0;
            <span class="keyword">end</span>
</pre><pre class="codeinput">        <span class="keyword">end</span>


    <span class="keyword">end</span>
</pre><pre class="codeinput"><span class="keyword">end</span>
</pre><h2>Part #2 , Checking the stability of intereset points using Laplacian of Gaussian (LOG)<a name="30"></a></h2><p>Choosing the neighbours and assembling corner points for all the scales</p><pre class="codeinput">[ro1, col1] = find(myLogic);
</pre><p>LOG of the image at a scale</p><pre class="codeinput">LOG{n} = (sigma^2)*sqrt(grad_xx.^2 + grad_yy.^2);
L_current = LOG{n};
</pre><p>Assembling all the potential interset points ( outputs of R threshold) in an array ( 3 X N)</p><pre class="codeinput">index_log = ro1 + (col1 - 1)*size(myLogic,1);
myPoints{n} = [ro1, col1, sigma*ones(numel(ro1),1)];
<span class="comment">% figure, imshow(uint8(myImage));hold on,plot(col1, ro1,'+');</span>
</pre><p>Maximizing the LOG and choosing the final HARRIS points. for the two extreme scales ( Top and Bottom) we only have one neighbor to compare with, but for rest of the scales compare LOG accross three scales.</p><p>Uppermost scale</p><pre class="codeinput"><span class="keyword">if</span> n == 2
</pre><pre class="codeinput">    L = LOG{1};
    myP = myPoints{1};
    index_check = myP(:,1) + (myP(:,2) - 1)*size(myLogic,1);
    diff_2 = L(index_check) - L_current(index_check);
</pre><p>Check for the maximum range</p><pre class="codeinput">    <span class="keyword">for</span> check = 1: size(myPoints{1},1)
        <span class="keyword">if</span> diff_2(check)&gt;0
            my_final = [my_final ; myP(check,1),myP(check,2),myP(check,3)];
        <span class="keyword">end</span>
    <span class="keyword">end</span>
</pre><p>All the middle range scales</p><pre class="codeinput"><span class="keyword">elseif</span> n&gt;2 &amp;&amp; n&lt;levels
</pre><pre class="codeinput">    L_center = LOG{n-1};
    L_up = LOG{n-2};
    myP = myPoints{n-1};
    index_check = myP(:,1) + (myP(:,2) - 1)*size(myLogic,1);
    diff_down = L_center(index_check) - L_current(index_check);
    diff_up = L_center(index_check) - L_up(index_check);
     <span class="keyword">for</span> check = 1: size(myP,1)
        <span class="keyword">if</span> diff_down(check)&gt;0 &amp;&amp; diff_up(check) &gt;0
            my_final = [my_final ; myP(check,1),myP(check,2),myP(check,3)];
        <span class="keyword">end</span>
    <span class="keyword">end</span>
</pre><p>Bottom scale</p><pre class="codeinput"><span class="keyword">elseif</span> n==levels
        L_up = LOG{n-1};
        myP = myPoints{10};
        index_check = myP(:,1) + (myP(:,2) - 1)*size(myLogic,1);
        diff_10 = L_current(index_check) - L_up(index_check);
    <span class="keyword">for</span> check = 1: size(myP,1)
         <span class="keyword">if</span> diff_10(check)&gt;0
            my_final = [my_final ; myP(check,1),myP(check,2),myP(check,3)];
        <span class="keyword">end</span>
    <span class="keyword">end</span>

<span class="keyword">end</span>
</pre><pre class="codeinput"><span class="keyword">end</span>
<span class="comment">% figure, imshow(uint8(imread(Img)));hold on;</span>
<span class="comment">% plot(my_final(:,2), my_final(:,1),'r*');</span>
</pre><p class="footer"><br><a href="http://www.mathworks.com/products/matlab/">Published with MATLAB® R2014a</a><br></p></div>

<div class="content"><h1></h1><!--introduction--><!--/introduction--><h2>Contents</h2><div><ul><li><a href="#1">Computation of HOMOGRAPHY using RANSAC</a></li><li><a href="#2">References:</a></li></ul></div><h2>Computation of HOMOGRAPHY using RANSAC<a name="1"></a></h2><h2>References:<a name="2"></a></h2><div><ol><li>Digital Video Processing (<a href="http://10.6.4.152/dvp/dvp.html">http://10.6.4.152/dvp/dvp.html</a>), Computer Sc., Indian Institute of Technology, Madras.</li><li>VLFEAT SIFT Tool box (<a href="http://www.vlfeat.org/overview/sift.html">http://www.vlfeat.org/overview/sift.html</a>)</li><li>RANSAC algorithm with example of finding homography : Edward Wiggin. MATLAB Central 2011.</li></ol></div><pre class="codeinput"><span class="keyword">function</span> H = RANSAC(sift_match_points_f1, sift_match_points_f2)
lmax = size(sift_match_points_f1, 2);
C = 0;

<span class="keyword">while</span> C &lt; 0.85*lmax
<span class="comment">% Picking the Random points (4 Nos)</span>
    ind = randIndex(lmax, 4);
    pts1 = sift_match_points_f1(:,ind);
    pts2 = sift_match_points_f2(:,ind);

<span class="comment">% Compute the HOMOGRAPHY</span>
    H = solveHomo(pts1,pts2);

<span class="comment">% Compute the residual for rest of points (lmax - 4)</span>
    remain_pts = H*[sift_match_points_f1;ones(1,lmax)];
    norm_remain_pts = [remain_pts(1,:)./remain_pts(3,:);<span class="keyword">...</span>
                      remain_pts(2,:)./remain_pts(3,:);ones(1,lmax)];

    distance_ransac  = norm_remain_pts - [sift_match_points_f2;<span class="keyword">...</span>
                                                         ones(1,lmax)];

    magnitude_distance = sqrt(sum(distance_ransac.^2,1));
    myRansacLogic = magnitude_distance &lt; 10;
    C = size(find(myRansacLogic),2);

<span class="keyword">end</span>
</pre><pre class="codeoutput error">Error using RANSAC (line 11)
Not enough input arguments.
</pre><p class="footer"><br><a href="http://www.mathworks.com/products/matlab/">Published with MATLAB® R2014a</a><br></p></div>


