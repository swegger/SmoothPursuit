function x = MultiSizeSensorModel(z,H,Wx,varargin)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'z')
addRequired(Parser,'H')
addRequired(Parser,'Wx')

parse(Parser,z,H,Wx,varargin{:})

z = Parser.Results.z;
H = Parser.Results.H;
Wx = Parser.Results.Wx;

%% Generate model sensor responses
eta = repmat(z,[size(H,1),1]).*repmat(Wx,[1,length(z)]).*randn(size(H,1),length(z));
x = H*z + eta;