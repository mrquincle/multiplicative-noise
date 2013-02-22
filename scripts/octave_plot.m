if (nargin < 1)
    printf('Usage: octave %s FILENAME (but now %d arguments)\n\n', program_name(), nargin); 
    exit
endif

arg_list = argv();
for i = 1:nargin
    printf('%d argument = %s\n', i, arg_list{i});
endfor 

file = arg_list{1}

a=load(file);
t=a(:,1);
x=a(:,2);

cut = length(t); % * 0.4
t=t(1:cut);
x=x(1:cut);

lt=log10(t);
lx=log10(x);

close all;
grid on;

%loglog(t,x);

plot(lt,lx);

Y = lx;
X = [lt ones(length(lt),1)];

% solving for m and c
alpha = inv(X'*X)*X'*Y

% constructing the straight line using the estimated slope and constant
yEst = alpha(1)*lt + alpha(2);

hold on;
plot(lt,yEst)
sleep(2)
input ("Hit key to exit")
