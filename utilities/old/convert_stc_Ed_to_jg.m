function convert_stc_Ed_to_jg(Trial,stcMode,nRev,varargin)
if isempty(stcMode),
    stcMode = 'manual1';
end

if ~isa(Trial,'MTASession'),
Trial = MTATrial.validate(Trial);
end
states = {'walk','rear','turn','pause','groom','sit','shake'};
keys   = {'w'   ,'r'   ,'n'   ,'p'    ,'m'    ,'s'  ,'k'};   

Trial.load('stc',stcMode);

while numel(varargin)>0,    
    switch varargin{1}
      case 'convert_keys'

        for s = 1:numel(states),
            Trial.stc.states{Trial.stc.gsi(states{s})}.key = keys{s};
        end

      case 'mutex_states'
        Stc = Trial.stc.copy;
        
        Trial.stc.states{Trial.stc.gsi('turn')}.data =  Stc{'n-r-m-s'}.data;
        Trial.stc.states{Trial.stc.gsi('walk')}.data =  Stc{'w-n-r'}.data;
        Trial.stc.states{Trial.stc.gsi('pause')}.data = Stc{'p-n-r-w'}.data;
        Trial.stc.states{Trial.stc.gsi('groom')}.data = Stc{'m-p-s-r-w'}.data;    
        Trial.stc.states{Trial.stc.gsi('sit')}.data =   Stc{'s-m-p-r-w'}.data;    
      
    end
    varargin(1) = [];
end


Trial.stc.updateMode(['hand_labeled_rev' num2str(nRev) '_Ed']);
Trial.stc.save(1);
