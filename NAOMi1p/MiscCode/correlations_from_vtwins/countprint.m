function out = countprint(x)
% prints \b x, returns string for printing this as an output if nargout ==
% 1

if x == 1
  out = num2str(x);
else
  out = [repmat('\b',1,length(num2str(x-1))) num2str(x)];
end
if nargout == 0
  fprintf(out)
end