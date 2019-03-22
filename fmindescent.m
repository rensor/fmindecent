classdef fmindecent < handle
  
  properties(Constant)
    name = 'fmindescent';
  end

  properties
    options = [];
  end
  
  properties(SetAccess = public, Hidden = true)
    
    % Inputs
    fun = [];
    x0 = [];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    nonlcon = [];
    nDV = [];
    
    % Merit function variables
    f0 = []; % Initial, unpeanalized objective function value
    aFac = []; % Scaling parameter for merit function
    lambda = []; % Lagrange variables for "Augmented Lagrange (AL) method"
    nG = []; % Total number of constraints
    Asum = []; % Sum all columns in A
    AeqSum = []; % Sum all columns in Aeq
    g = []; % all constraints at latest evaluation. Stores as g <=0 in the following orderfields
            % [gnl;gnleq;-gnleq;glin;glineq;-glineq;-x+this.lb;x-this.ub];
            
    % Second order Quasi-Newton Methods
    Hinv = []; % Inverse of hessian DFP method
    H = []; % Hessian BFGS method
    iH = 0; % Hessian update counter
    
    % Iteration history
    history = struct();
    
    % Switch
    initialized = false;
    
    
  end
  
  methods
    
    % Construct
    function this = fmindecent(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,varargin)
      
      % initialize data structures
      if nargin < 1 || ~isa(fun,'function_handle')
        error([this.name,': fun is required to be a function handle'])
      else
        this.fun = fun;
      end
      
      if nargin < 2 || isempty(x0)
        this.x0 = [];
      else
        if isnumeric(x0)
          this.x0 = x0(:);
        else
          error([this.name,' x0 is required to be of type numeric'])
        end
      end
      
      if nargin < 3 || isempty(A)
        this.A = [];
      else
        if isnumeric(A)
          this.A = A;
          this.Asum = sum(A,2);
        else
          error([this.name,' A is required to be of type numeric'])
        end
      end
      
      if nargin < 4 || isempty(b)
        this.b = [];
      else
        if isnumeric(b)
          this.b = b(:);
        else
          error([this.name,' b is required to be of type numeric'])
        end
      end
      
      if size(this.A,1) ~= size(this.b,1)
        error([this.name,' A and b must contain equal number of rows'])
      end
      
      if nargin < 5 || isempty(Aeq)
        this.Aeq = [];
      else
        if isnumeric(Aeq)
          this.Aeq = Aeq;
          this.AeqSum = sum(Aeq,2);
        else
          error([this.name,' Aeq is required to be of type numeric'])
        end
      end
      
      if nargin < 6 || isempty(beq)
        this.beq = [];
      else
        if isnumeric(beq)
          this.beq = beq(:);
        else
          error([this.name,' beq is required to be of type numeric'])
        end
      end
      
      if size(this.Aeq,1) ~= size(this.beq,1)
        error([this.name,' Aeq and beq must contain equal number of rows'])
      end
          
      if nargin < 7 || isempty(lb)
        this.lb = [];
      else
        if isnumeric(lb)
          this.lb = lb(:);
        else
          error([this.name,' lb (lower bound) is required to be of type numeric'])
        end
      end
      
      if nargin < 8 || isempty(ub)
        this.ub = [];
      else
        if isnumeric(ub)
          this.ub = ub(:);
        else
          error([this.name,' ub (upper bound) is required to be of type numeric'])
        end
      end
      
      if nargin < 9 || isempty(nonlcon)
        this.nonlcon = [];
      elseif ~isempty(nonlcon)
        if ~isa(nonlcon,'function_handle')
          error([this.name,' nonlcon is required to be a function handle'])
        else
          this.nonlcon = nonlcon;
        end
      end
      
      % check that all sizes match
      if ~isempty(this.x0)
        this.nDV = numel(this.x0);
      end
      
      if ~isempty(this.lb) && ~isempty(this.ub)
        if numel(this.lb) ~= numel(this.ub)
          error([this.name,' lb and ub not equal dimensions'])
        end
      end
      
      if ~isempty(this.lb) && ~isempty(this.nDV)
        if numel(this.lb) ~= this.nDV
          error([this.name,' x0 and lb not equal dimensions'])
        end
      end
      
      if ~isempty(this.ub) && ~isempty(this.nDV)
        if numel(this.ub) ~= this.nDV
          error([this.name,' x0 and ub not equal dimensions'])
        end 
      end
      
      if ~isempty(this.A) && ~isempty(this.nDV)
        if size(this.A,2) ~= this.nDV
          error([this.name,' Columns of A(',num2str(size(this.A,2)),') does not match number of design variables(',num2str(this.nDV),')'])
        end
      elseif ~isempty(this.A) && isempty(this.nDV)
        this.nDV = size(this.A,2);
      end
      
      if ~isempty(this.Aeq) && ~isempty(this.nDV)
        if size(this.Aeq,2) ~= this.nDV
          error([this.name,' Columns of Aeq(',num2str(size(this.A,2)),') does not match number of design variables(',num2str(this.nDV),')'])
        end
      elseif ~isempty(this.Aeq) && isempty(this.nDV)
        this.nDV = size(this.Aeq,2);
      end
      
      % initialize options structure
      this.options = fmindecent.setOptions(varargin);
      
      % We made it this far
      this.initialized = true;
    end
    
    % Main function
    function [x,fval,exitflag,output] = solve(this)
      % Assume the code failed
      exitflag = -1;
      
      if strcmpi(this.options.Display,'iter')
        fprintf('*********************************************************************************************************')
        fprintf('\n \t \t \t \t \t \t fmindecent optimizer')
        fprintf('\n*********************************************************************************************************\n')
        fprintf('\t %10s \t\t %10s \t \t %10s \t \t   %10s \t \t   %10s \n','f(x)','Max inf', 'Norm dx', 'nFeval','IterNo');
      end
        
      % Allocate iteration history array
      % Store function values, maximum infeasibility from non-linear
      % constraints, norm of design variable change
      this.history.f = zeros(this.options.MaxIterations,1);
      this.history.xnorm = zeros(this.options.MaxIterations,1);
      if ~isempty(this.nonlcon)
        this.history.maxInf = zeros(this.options.MaxIterations,1);
      else
        % Assumed allocated to empty when not active
        this.history.maxInf = [];
      end
      
      % Ensure that we have a column vector (nDV,1)
      x = this.x0(:);

      % evaluate non-linear constraints
      if ~isempty(this.nonlcon)
        [gnl, gnleq] = this.nonlcon(x);
      else
        gnl = [];
        gnleq = [];
      end
      % evaluate linear in-equality constraints
      if isempty(this.A)
        glin = [];
      else
        glin = this.A*x - this.b;
      end
      
      % evaluate linear equality constraints
      if isempty(this.Aeq)
        glineq = [];
      else
        glineq = this.Aeq*x - this.beq;
      end
      % Assemble all constraints as leq<=0 type
      this.g = [gnl;gnleq;-gnleq;glin;glineq;-glineq;-x+this.lb;x-this.ub];
      % Count total number of constraints
      this.nG = size(this.g,1);
      % Determine initial infeasibility
      maxInf = max([this.g;0]);
      
      if strcmpi(this.options.ConstraintMethod,'AL')  
        % Initialize lagrange multipliers to zeros
        this.lambda = zeros(this.nG,1);
        % Update lagrange multipliers based on the initial infeasibilities
        this.updateLambda();
      end
      
      % evaluate objective function at initial point
      this.f0  = this.fun(x);
      
      % Get scaling factor for merit function
      this.aFac = max([abs(this.f0),1]);
      
      % Evaluate merit function
      fmerit = this.getMeritObj(x,false);
      
      % Store initial objective function value
      fOld = fmerit;
      
      % Create empty variable for "old" gradients
      dfmerit = [];
      dc = [];
      alpha = [];
      
      % Set counters and switches
      nFeval = 1;
      iterNo = 0;
      optimize = true;
      
      if strcmpi(this.options.Display,'iter')
        fprintf('\t %6.4e \t \t %6.4e \t \t %6.4e \t \t %10i \t \t %10i \n' ,this.f0, maxInf, 0, nFeval ,iterNo);
      end
        
      % Main loop
      while optimize
        
        % update iteration counter
        iterNo = iterNo + 1;
        
        if strcmpi(this.options.ConstraintMethod,'AL')  
        % update lagrange multipliers
          this.updateLambda()
        end
        
        % store previous gradients
        dfmeritOld = dfmerit;
        
        % evaluate gradients
        [~,dfmerit,~,~,dg] = this.getMeritObj(x,true);
        
        % Get decent direction based on user defined algorithm
        [dc,ndc] = this.getDecentDirection(dfmerit,dfmeritOld,dc,alpha,dg);
        
        % call linear search method
        if ndc > 0
          [xNew,alpha,fmerit,fval,nF,exitflag] = this.lineSearch(dc,x);
        else
          xNew = x;
          alpha = 0;
          fmerit = fOld;
          nF = 0;
        end
      
        % Update function evaluation counter
        nFeval = nFeval + nF;
        
        maxInf = max([this.g;0]);
        
        optimalityNorm = sqrt((fOld-fmerit)^2);
        % check for convergence
        if ( (optimalityNorm <= this.options.OptimalityTolerance) || (alpha <=this.options.StepTolerance)) || (iterNo >= this.options.MaxIterations) || (nFeval >= this.options.MaxFunctionEvaluations)
          optimize = false;
          exitflag = 1;
        end
        
        % Update design variables
        x = xNew;
        
        % Update "old" design
        fOld = fmerit;
        
        % Store iteration history
        this.history.f(iterNo) = fval;
        this.history.xnorm(iterNo) = alpha;
        this.history.maxInf(iterNo) = maxInf;
        this.history.nIter = iterNo;
        this.history.nFeval = nFeval;
        
        if strcmpi(this.options.Display,'iter')
            fprintf('\t %6.4e \t \t %6.4e \t \t %6.4e \t \t %10i \t \t %10i \n' ,fval, maxInf, alpha, nFeval ,iterNo);
        end
        
      end % Main loop
      
      this.history.f(iterNo+1:end)=[];
      this.history.xnorm(iterNo+1:end)=[];
      this.history.maxInf(iterNo+1:end)=[];
      output.history = this.history;
      if strcmpi(this.options.ConstraintMethod,'AL')
        output.lambda = this.lambda./this.aFac;
      end
      
    end % Solve function
    
    function postprocess(this)
      % Save current "default" window style
      defaultWindowStyle=get(0,'DefaultFigureWindowStyle');
      % Set new window style to docked
      set(0,'DefaultFigureWindowStyle','docked')
      
      % Make iteration vector
      ivec = 1:this.history.nIter;
      
      f1=figure();
      plot(ivec,this.history.f)
      title('Objective')
      xlabel('Iteration Number')
      ylabel('Objective value')
      
      figure();
      plot(ivec,this.history.xnorm)
      title('Design change norm')
      xlabel('Iteration Number')
      yl=ylabel('Norm dx');
      set(yl,'Interpreter','none')
      
      figure();
      plot(ivec,this.history.maxInf)
      title('Maximum infeasibility')
      xlabel('Iteration Number')
      ylabel('-')
      
      % Jump back to figure 1
      figure(f1)
      % Restore default window style
      set(0,'DefaultFigureWindowStyle',defaultWindowStyle)
    end
    
  end % methods
  
  
  methods (Hidden = true)
    
    function [fmerit,dfmerit,fval,df,dgActive] = getMeritObj(this,x,doDSA)
      fmerit = [];
      fval = [];
      dfmerit = [];
      df = [];
      dgActive = [];
      if nargin < 3 || isempty(doDSA)
        doDSA = false;
      end
      
      if ~doDSA
        fval = this.fun(x);
        % evaluate non-linear constraints
        if ~isempty(this.nonlcon)
          [gnl, gnleq] = this.nonlcon(x);
        else
          gnl = [];
          gnleq = [];
        end
        
        % evaluate linear in-equality constraints
        if isempty(this.A)
          glin = [];
        else
          glin = this.A*x - this.b;
        end
        % evaluate linear equality constraints
        if isempty(this.Aeq)
          glineq = [];
        else
          glineq = this.Aeq*x - this.beq ;
        end
        % Assemble all constraints as leq<=0 type
        this.g = [gnl;gnleq;-gnleq;glin;glineq;-glineq;-x+this.lb;x-this.ub];
        % determine infeasibility
        y = max(this.g ,0);
        switch this.options.ConstraintMethod
          case 'Merit'
            fmerit = fval + this.aFac*sum(y*this.options.InfeasibilityPenalization+0.5*y.^2);
          case 'AL'
            fmerit = fval + this.aFac*sum(y.*this.lambda+this.options.InfeasibilityPenalization*0.5*y.^2);
        end
      else
          [~,df] = this.fun(x);
          
          if ~isempty(this.nonlcon)
            [~,~,dgnl,dgneq] = this.nonlcon(x);
          else
            dgnl = [];
            dgneq = [];
          end

          dgreal = [dgnl';dgneq';-dgneq'; this.A; this.Aeq;-this.Aeq; -eye(this.nDV); eye(this.nDV)];
          % Apply active set strategy
          Active = this.g>=0;
          dgActive = dgreal(Active,:);
          switch this.options.ConstraintMethod
            case 'Merit'
              % Sum sensitivites from all constraints together for each design variable
              temp = sum(dgreal(Active),1)';
              dg = this.aFac.*(this.options.InfeasibilityPenalization+temp);
            case 'AL'
              dg = sum(this.aFac.*(this.lambda(Active)+this.options.InfeasibilityPenalization.*dgreal(Active,:)),1)';
          end
          
          
          % Add sensitivites from objective and constraints
          dfmerit = df+dg;
      end
    end
    
    function updateLambda(this)
      this.lambda = max(this.lambda + this.options.InfeasibilityPenalization*this.g,0);
    end
    
    function [dc,ndc] = getDecentDirection(this,df,dfm1,dcm1,alpha,dg)
      switch this.options.Algorithm
        case 'CG'
          % Conjugate gradient update: Fletcher-Reeves Method
          if ~isempty(dcm1) 
            beta = df'*df/(dfm1'*dfm1); % from second iteration, use CG
            dc = -df + beta*dcm1;
          else
            dc = -df; % first iteration, use stepest decent
          end
        case 'DFP' 
          % DFP inverse hessian update (Quasi Newton)
          if ~isempty(dcm1) 
            s = alpha.*dcm1;
            y = df-dfm1;
            z = this.Hinv*y;
            C = -z*z'/(dot(y,z));
            B = s*s'/(dot(s,y));
            if this.iH >= this.options.HessianRest
              this.Hinv = eye(this.nDV);
              this.iH = 0;
            end
            this.iH = this.iH + 1;
            this.Hinv = this.Hinv + B + C;
            dc = -this.Hinv*df;
          else
            this.Hinv = eye(this.nDV);
            dc = -df;
          end
        case 'BFGS'
          % BFGS hessian update (Quasi Newton)
          if ~isempty(dcm1) 
            
            if this.iH >= this.options.HessianRest
              this.H = eye(this.nDV);
              this.iH = 0;
            end
            
            s = alpha.*dcm1;
            y = df-dfm1;
            sHs = s'*this.H*s;
            sy = s'*y;
            if sy >=0.2*sHs
              theta = 1;
              r = y;
            else
              theta = 0.8*sHs/(sHs-sy);
              r = theta*y+(1-theta)*this.H*s;
            end
            this.iH = this.iH + 1;
            this.H = this.H - this.H*s*s'*this.H/(sHs)+r*r'/(s'*r);
            dc = -this.H/df';
          else
            this.H = eye(this.nDV);
            dc = -df;
          end
          
        case 'MFD'
          
          % Get number of active constraints
          ng = size(dg,1);
          if ng > 0
            % Push off parameter
            theta = 0.8;
            % Define variable type
            ndv = numel(df)+1;
            vartype = char(ndv,1);
            vartype(1:ndv) = 'C';
            % Define upper and lower bounds
            mfdlb = -ones(ndv,1);
            mfdlb(end) = 0; % beta variable
            mfdub = ones(ndv,1);
            mfdub(end) = 1000; % beta variable
            
            % Specify obj gradient for MFD problem
            dBeta = zeros(ndv,1);
            dBeta(end) = 1; % beta variable
            
            % Define constraints
            
            mfdA = zeros(ng+1,ndv);
            mfdA(1,1:end-1) = df;
            mfdA(2:end,1:end-1) = dg;
            mfdA(1,end) = -1; % beta variable
            mfdA(2:end,end) = -1*theta; % beta variable
            mfdB = zeros(ng+1,1);
            
            nleq = ng+1;
            ctype = char(nleq,1);
            ctype(1:nleq) = 'U';
            
            [temp, ~, exitflag] = glpk (dBeta, mfdA, mfdB, mfdlb, mfdub, ctype, vartype);
            dc = temp(1:end-1);
          else
            dc = -df;
          end
          
        otherwise
          dc = -df;
      end
      ndc = norm(dc);
      if ndc > 0
        dc = dc./norm(dc); % Normalize 
      end
      
    end
    
    function [xNew, alpha,fmerit,fval,nF,exitflag] = lineSearch(this,df,x)
      switch this.options.LineSearch
        case 'golden'
          [xNew, alpha,fmerit,fval,nF,exitflag] = this.goldenSectionSearch(df,x);
      end
        
    end
    
    function [xNew,alpha,fmerit,fval,nF,exitflag] = goldenSectionSearch(this,dc,x)
      
      % Define "golden" constants
      phi = (sqrt(5)+1)/2;
      invPhi = 1/phi;
      invPhi2 = 1/phi^2;
      
      % Set function evaluation counter
      nF = 0;
      
      % Initial step length along decent direction
      delta = 0.001; 
      % Lower limit on alpha (step lenght)
      alphaL = 0; 
      
      % Initial bracketing
      yl = this.getMeritObj(x);
      alphaU = 0; % Initialize upper limit on step length
      alphaUm1 = 0; % Initialize upper limit minus 1
      for ii = 0:this.options.MaxFunctionEvaluations-1
        alphaUm2 = alphaUm1;
        alphaUm1 = alphaU;
        alphaU = alphaU + delta * phi^ii;
        xu = x + dc*alphaU;
        nF = nF + 1;
        yu = this.getMeritObj(xu);
        if yu > yl
          alphaL=alphaUm2;
          break
        end
        yl = yu;
      end
        
      % Interval reduction
      h = alphaU-alphaL;
      
      if h <= this.options.StepTolerance
        alpha = (alphaL+alphaU)/2;
        xNew = x + dc*alpha;
        nF = nF + 1;
        [fmerit,~,fval] = this.getMeritObj(xNew);
        exitflag = 1;
        return
      end
      
      % required steps to reach tolerance
      n = ceil(log(this.options.StepTolerance/h)/log(invPhi));
      
      % Define inner interval
      c = alphaL + invPhi2*h;
      d = alphaL + invPhi*h;
      % Evaluate point c
      xc = x + dc*c;
      nF = nF + 1;
      yc = this.getMeritObj(xc);
      % Evaluate point d
      xd = x + dc*d;
      nF = nF + 1;
      yd = this.getMeritObj(xd);
      
      for ii = 1:n
        if yc < yd
          alphaU = d;
          d = c;
          yd = yc;
          h = invPhi*h;
          c = alphaL + invPhi2*h;
          % Evaluate point c
          xc = x + dc*c;
          nF = nF + 1;
          yc = this.getMeritObj(xc);
        else
          alphaL = c;
          c = d;
          yc = yd;
          h = invPhi*h; 
          d = alphaL + invPhi*h;
          % Evaluate point d
          xd = x + dc*d;
          nF = nF + 1;
          yd = this.getMeritObj(xd);
        end
      end
      
      if yc < yd
        c = alphaL;
        alpha = (c+d)/2;
      else
        d = alphaU;
        alpha = (c+d)/2;
      end
      
      % Evaluate final point
      xNew = x + dc*alpha;
      nF = nF + 1;
      [fmerit,~,fval] = this.getMeritObj(xNew);
      exitflag = 1;
    end
    
  end
  
  methods(Static = true, Hidden = true)
    
    % options initialization
    function options = setOptions(input)
      % Here you can add new options if needed
      p = inputParser;
      p.CaseSensitive = false;
      % Helper functions for input parser
      checkEmpetyOrChar = @(x) (isempty(x) || ischar(x));
      checkEmptyOrNumericPositive = @(x) (isempty(x) || (isnumeric(x) && all(x > 0)));
      
      % Set parameters
      p.addParameter('Algorithm','CG',  @(x) checkEmpetyOrChar(x));
      p.addParameter('ConstraintMethod','AL',  @(x) checkEmpetyOrChar(x));
      p.addParameter('LineSearch','golden',  @(x) checkEmpetyOrChar(x));
      p.addParameter('Display','off',  @(x) checkEmpetyOrChar(x));
      p.addParameter('MaxFunctionEvaluations',1000,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('MaxIterations',1000,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('InfeasibilityPenalization',100,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('OptimalityTolerance',1e-5,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('StepTolerance',1e-5,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('HessianRest',15,  @(x) checkEmptyOrNumericPositive(x));
      
      % pars input
      if nargin < 1 || isempty(input)
        parse(p);
      else
        parse(p,input{:});
      end
      
      % Output results to options structure
      options = p.Results;
      
    end
  end
end