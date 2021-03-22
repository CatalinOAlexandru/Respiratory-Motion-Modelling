% Displaying all images and make an animation of MR images for task 6.1

for k = 1:100
  disp(k);
  
  % Runs only through first 100 images as its a little slow.
  % Problably because it has to handle many windows open and closes.
  
  % Keep the next line otherwise matlab cant save the figure.
  f = figure;
  dispNiiSlice(images(k),'z',1);
  drawnow;
  

  % Stole this code from google.  
  frame = getframe(f); 
  im = frame2im(frame); 
  [imind,cm] = rgb2ind(im,256); 
  % Write to the GIF File 
  if k == 1 
      imwrite(imind,cm,'test.gif','gif', 'Loopcount',inf); 
  else 
      imwrite(imind,cm,'test.gif','gif','WriteMode','append'); 
  end 
  clf
  close % keep this 'close' otherwise it wont close any figure window and you
  % will end up with N number of overlapping windows.
end


                                   
                                   
                                   
                                   